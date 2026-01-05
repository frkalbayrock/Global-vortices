module Time_Evolution

include("Parameters.jl")
include("Energy.jl")
include("IndexMap.jl")
include("Auxiliary.jl")
include("Constraints.jl")
include("FindVortex.jl")
using Distributed
using .Parameters 
using .Energy
using .IndexMap
using .Auxiliary_Routines
using .Constraints_Conserveds
using .FindVortex
using DelimitedFiles





# 4-index Notation (Leap-frog)
export time_evolve!
function time_evolve!(ϕ_gl,ψ_gl,ZED_gl,meanSqrRenorm,zPE)

    #--Snapshotting interval
    if round(nt/nsnaps,RoundDown) == 0
        snapInterval = 1
        println("---Caution: Number of time steps is smaller than snaps!---")
        println("Instead of ",nsnaps," snapshots, ",nt," snaps will be taken.","---")
    else
        snapInterval = round(Int,nt/nsnaps,RoundDown)
        if nt/nsnaps != snapInterval
            actualSnaps=round(Int,nt/snapInterval,RoundDown)
            println("---Caution: Number of snapshots are not multiple of time steps.")
            println("Instead of ",nsnaps," snapshots, ",actualSnaps," snaps will be taken.","---")
        end
    end


    #--Open data files
    ϕdataIO = open("data/phi.dat","w") 
    ψdataIO = open("data/psi.dat","w")
    energyIO = open("data/energy.dat","w")
    ZedIO = open("data/energies/ZED.dat","w")
    vortexIO = open("data/vortices.dat","w")




    #-Time Evolution
    for t=1:nt

        #--Move a time step
        @everywhere workers() half_step!(ϕ,ψ,Z,dϕdt,dψdt,dZdt,1)
        @everywhere workers() update_Paddings!(ϕ,ψ,Z)
        @everywhere workers() begin
            leap_forward!(ϕ,ψ,Z,dϕdt,dψdt,dZdt,($meanSqrRenorm))
            half_step!(ϕ,ψ,Z,dϕdt,dψdt,dZdt,2)
        #Shift next time values to present time for the next step
            updateForNextStep(dϕdt,dψdt,dZdt)
        end


        #--Take a snap
        if mod(t,snapInterval) == 0

            #-Calculate energy
            totalE = energy!(ZED_gl,meanSqrRenorm,zPE)
            

            #-Update global fields for recording
            for i=2:nprocs()
                lx_p, rx_p, ly_p, ry_p = chunker(i) #_p: physical
                ϕ_gl[lx_p:rx_p,ly_p:ry_p] .= (@fetchfrom i Main.ϕ[lx_l:rx_l,
                                                                  ly_l:ry_l])::Array{ComplexF64, 2}
                ψ_gl[lx_p:rx_p,ly_p:ry_p] .= (@fetchfrom i Main.ψ[lx_l:rx_l,
                                                                  ly_l:ry_l])::Array{Float64, 2}
            end


            #----Check for Vortices----#
            vortex_pos, anti_vortex_pos = vortex_finder(ϕ_gl,t)


            #Record vortices
            if length(vortex_pos) != length(anti_vortex_pos)
                error("Number of vortices doesn't match anti-vortices!")
            elseif (length(vortex_pos)==0 && length(anti_vortex_pos)==0)
                println(vortexIO,"[]")   #vortex 
                println(vortexIO,"[]")   #anti-vortex
            else
                println(vortexIO,vortex_pos)         #vortex
                println(vortexIO,anti_vortex_pos)    #anti-vortex
            end


            #Record field and energy data
            writedlm(ϕdataIO, ϕ_gl[:,:])
            writedlm(ψdataIO, ψ_gl[:,:])
            writedlm(energyIO,totalE)
            writedlm(ZedIO,ZED_gl)


            #!
            #Check constraints
            # Z_gl_f = mapZTo2Index(Z_gl[:,:,:,:])
            # dZdt_gl_f = mapZTo2Index(dZdt_gl[:,:,:,:])
            # @views constraints_checker(Z_gl_f,dZdt_gl_f)
            # @views conserved_checker(Z_gl_f,dZdt_gl_f)

        end 
    end

    close(ϕdataIO)
    close(ψdataIO)
    close(energyIO)
    close(ZedIO)
    close(vortexIO)
    
end






#Half-step (vectorized)
@everywhere workers() function half_step!(ϕ,ψ,Z,dϕdt,dψdt,dZdt,t)
    @views begin
    @. ϕ[padd+1:Nx_loc+padd,padd+1:Ny_loc+padd]     += @views (dt_half) * dϕdt[:,:,t]
    @. ψ[padd+1:Nx_loc+padd,padd+1:Ny_loc+padd]     += @views (dt_half) * dψdt[:,:,t]
    @. Z[padd+1:Nx_loc+padd,:,padd+1:Ny_loc+padd,:] +=  (dt_half) * @views dZdt[:,:,:,:,t]
    end
return nothing
end




#Leap-forward (vectorized)
@everywhere workers() function leap_forward!(ϕ,ψ,Z,dϕdt,dψdt,dZdt,meanSqrRenorm)

    fluxes_ϕ_ψ!(ϕ,ψ,Z,ϕ_flux,ψ_flux,meanSqr_Rho,meanSqrRenorm)
    flux_Z!(ϕ,ψ,Z,Z_flux)
    @. @views dϕdt[:, :, 2] = dϕdt[:, :, 1] + dt*( ϕ_flux )
    @. @views dψdt[:, :, 2] = dψdt[:, :, 1] + dt*( ψ_flux )
    @. @views dZdt[:,:,:,:,2] = dZdt[:,:,:,:,1] + dt*( Z_flux )

return nothing
end




#Updating the time coordinates for the next step
@everywhere workers() function updateForNextStep(dϕdt,dψdt,dZdt)

    #Shift time coordinates
    @views begin
        dϕdt[:,:,1] .= dϕdt[:,:,2]
        dψdt[:,:,1] .= dψdt[:,:,2]
        dZdt[:,:,:,:,1] .= dZdt[:,:,:,:,2]
    end

return nothing
end




#Fluxes for ϕ and ψ
@everywhere workers() @inline function fluxes_ϕ_ψ!(ϕ,ψ,Z,ϕ_flux,ψ_flux,meanSqr_Rho,meanSqrRenorm)

    #Save ranges that are gonna be used for clarity 
    j = padd+1 : Nx_loc+padd     
    k = padd+1 : Ny_loc+padd     

    @views begin
        #2-point function
        meanSqr_Rho .= @views dropdims(sum(abs2, Z[j, :, k, :]; dims=(2,4)), dims=(2,4))

        @. ϕ_flux = ( 
                      (ϕ[j+1, k] - 2ϕ[j, k] +  ϕ[j-1, k])/dx2 
                    + (ϕ[j, k+1] - 2ϕ[j, k] +  ϕ[j, k-1])/dy2
                    + ( m_ϕ^2 - α/2 *(meanSqr_Rho - meanSqrRenorm)/(dx*dy) - λ * abs2(ϕ[j,k]) ) * ϕ[j,k] 
                    )

        @. ψ_flux = ( 
                      (ψ[j+1, k] - 2ψ[j, k] +  ψ[j-1, k])/dx2 
                    + (ψ[j, k+1] - 2ψ[j, k] +  ψ[j, k-1])/dy2   
                    - ( m_ψ^2 + β *(meanSqr_Rho - meanSqrRenorm)/(dx*dy) ) * ψ[j,k] 
                    )
    end

end



#Flux for Z
@everywhere workers() @inline function flux_Z!(ϕ,ψ,Z,Z_flux)

    #Save ranges that are gonna be used for clarity 
    j = padd+1 : Nx_loc+padd     
    k = padd+1 : Ny_loc+padd     

    @views begin
        # ϕ[j,k] ψ[j,k]'s are 2D arrays so we expand them to 
        # 4D for vectorized calculation with Z: (j,1,k,1)
        M = reshape(@.(m_ρ^2 + α * abs2(ϕ[j, k]) + β * (ψ[j, k]^2)), 
                                                Nx_loc, 1, Ny_loc, 1)        

        # final vectorized update (broadcasting is explicit and correct)
        @. Z_flux = ( 
                    ( ( Z[j+1, :, k, :] - 2 * Z[j, :, k, :] +  Z[j-1, :, k, :]) / dx2 
                     +( Z[j, :, k+1, :] - 2 * Z[j, :, k, :] +  Z[j, :, k-1, :]) / dy2 ) 
                     - M * Z[j, :, k, :]
                    )
    end

end




end #module
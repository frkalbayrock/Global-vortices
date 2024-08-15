module Time_Evolution

include("Parameters.jl")
using .Parameters 
include("IndexMap.jl")
using .IndexMap
include("Auxiliary.jl")
using .Auxiliary_Routines
include("Constraints.jl")
using .Constraints_Conserveds
using DelimitedFiles
using OffsetArrays
using Profile
using PProf
using JLD2



export time_evolve!
#We use Leap-Frog: Position Verlet (LFPV) method to time evolve the system.
#Position: refers to fields (i.e. ψ(t)) , Velocity: refers to time derivatives (i.e. dψdt(t))
#LFPV first shifts the "position" to a half-integer time step with 
#half Euler step (integrating to t+1/2 using velocity at time t - forward finite diff.)
#Next, the "velocity" is leaped forward a full time step with 
#full Euler step using the midpoint(t+1/2) "position" calculated in the previous step.
#Then evolve the "position" back to an integer time-step using half-Euler now using "velocity" at the t+1 step.
function time_evolve!(ϕ,ψ,Z,dϕdt,dψdt,dZdt,meanSqrRenorm,zPE)

    #Flatten ϕ, ψ, dϕdt, dψdt #!no offset for these#!
    ϕ_s = im*zeros(N^2,2)
    ψ_s = zeros(N^2,2)
    dϕdt_s = im*zeros(N^2,2)
    dψdt_s = zeros(N^2,2)
    ϕ_s[:,1] =  @views flattenDimension(ϕ[:,:,0])   #!for these .= uses more allocations!
    ψ_s[:,1] =  @views flattenDimension(ψ[:,:,0])  
    dϕdt_s[:,1] =  @views flattenDimension(dϕdt[:,:,0])
    dψdt_s[:,1] =  @views flattenDimension(dψdt[:,:,0])

    #Snapshotting interval
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

    ϕdataIO = open("data/phi.dat","w") 
    ψdataIO = open("data/psi.dat","w")
    energyIO = open("data/energy.dat","w")
    ZedIO = open("data/energies/ZED.dat","w")

    # meanRhoSqrIO = open("data/meanRhoSqr.dat","w")
    # ϕJLD = jldopen("data/phi.jld2","a")#!
    # ψJLD = jldopen("data/psi.jld2","a")#!
    # count=0

    for t=1:nt
        
        #Move a time step
        half_step!(ϕ_s,ψ_s,Z,dϕdt_s,dψdt_s,dZdt,1)
        leap_forward!(ϕ_s,ψ_s,Z,dϕdt_s,dψdt_s,dZdt,meanSqrRenorm)
        half_step!(ϕ_s,ψ_s,Z,dϕdt_s,dψdt_s,dZdt,2)
        # Shift next time values to present time for the next step
        updateForNextStep!(ϕ_s,ψ_s,Z,dϕdt_s,dψdt_s,dZdt)

        # Take a snap
        if mod(t,snapInterval) == 0
            #!! UPDATE!!!
            #!these mappings are unnecessary just adjust your routines to work with single-index notation
            @views ϕ[:,:,0] .= ravelDimension(ϕ_s[:,1])    #! i just realized; if i'll just evolve single index ones i don't really need 
            @views ψ[:,:,0] .= ravelDimension(ψ_s[:,1])    #! time coordinates of original ϕ and ψ; just reduce them once the testing is over.
            @views dϕdt[:,:,0] .= ravelDimension(dϕdt_s[:,1])
            @views dψdt[:,:,0] .= ravelDimension(dψdt_s[:,1])

            totalE,ZED = energy(ϕ,ψ,Z,dϕdt,dψdt,dZdt,meanSqrRenorm,zPE)
            ZED_t = ravelDimension(ZED) #!unnecessary temp array.

            writedlm(ϕdataIO, ϕ[:,:,0])
            writedlm(ψdataIO, ψ[:,:,0])
            writedlm(energyIO,totalE)
            writedlm(ZedIO,ZED_t)
        

            #!
            # @views recordSnap(ϕ[:,:,0],ψ[:,:,0],totalE,ZED_t,
            #                     ϕdataIO,ψdataIO,energyIO,ZedIO)
    
            #!
            #check constraints
            # @views constraints_checker(Z[:,:,1],dZdt[:,:,1])
            # @views conserved_checker(Z[:,:,1],dZdt[:,:,1])
        end

    end
    # close(ϕJLD)#!
    # close(ψJLD)#!
end



#Half-step iterate "position" using present time values (t=0)
#Half-step iterate "position" using next time value of "velocity" (t=1)
function half_step!(ϕ,ψ,Z,dϕdt,dψdt,dZdt,t)
    # For t=0, fields are calculated at time 1/2 & for t=1 they are calculated fully at t=1.
    N2 = size(ϕ,1)  #! these seem to slightly improve allocations
    dt_half =dt/2   #!not sure yet if i wanna keep them. harder to read.
    # count =0#!
    Threads.@threads for J=1:N2
    # for J=1:N2
        # count +=1#!
        # println(count)#!
        ϕ[J,2] = ϕ[J,t] + dt_half*( dϕdt[J,t] )
        ψ[J,2] = ψ[J,t] + dt_half*( dψdt[J,t] )
        for K=1:N2
            Z[J,K,2] = Z[J,K,t] + dt_half*( dZdt[J,K,t] )
        end
    end

end


#Full-step iterate "velocity" using half-integer time value of "position" (t=1/2)
function leap_forward!(ϕ,ψ,Z,dϕdt,dψdt,dZdt,meanSqrRenorm)
    
    Threads.@threads for J=1:N^2
    # for J=1:N^2
        #PBC
        nnl_x, nnr_x, nnl_y, nnr_y = pbc1D(J)
        #Calculate fluxes for ϕ and ψ
        @views ϕ_flux , ψ_flux = fluxes_ϕ_ψ(ϕ[:,2],ψ[:,2],Z[:,:,2],meanSqrRenorm,J,nnl_x,nnr_x,nnl_y,nnr_y)
        dϕdt[J,2] = dϕdt[J,1] + dt*( ϕ_flux )
        dψdt[J,2] = dψdt[J,1] + dt*( ψ_flux )
        for K=1:N^2
            Z_flux = flux_Z(ϕ,ψ,Z,J,K,nnl_x,nnr_x,nnl_y,nnr_y)
            dZdt[J,K,2] = dZdt[J,K,1] + dt*( Z_flux )
        end
    end

end


function fluxes_ϕ_ψ(ϕ,ψ,Z,meanSqrRenorm,J,nnl_x,nnr_x,nnl_y,nnr_y)
    #2-point function
    meanSqr_Rho = @views sum(abs2, Z[J,:])

    ϕ_flux = ( (ϕ[nnr_x] - 2ϕ[J] + ϕ[nnl_x])/dx^2 + (ϕ[nnr_y] - 2ϕ[J] + ϕ[nnl_y])/dy^2 
            + ( m_ϕ^2 -α/2 *(meanSqr_Rho - meanSqrRenorm)/(dx*dy) -λ* abs2(ϕ[J]) ) * ϕ[J] )

    ψ_flux = ( (ψ[nnr_x] - 2ψ[J] + ψ[nnl_x])/dx^2 + (ψ[nnr_y] - 2ψ[J] + ψ[nnl_y])/dy^2
            - ( m_ψ^2 +β *(meanSqr_Rho - meanSqrRenorm)/(dx*dy) ) * ψ[J] )

return ϕ_flux, ψ_flux
end

function flux_Z(ϕ,ψ,Z,J,K,nnl_x,nnr_x,nnl_y,nnr_y)
    Z_flux = ( (Z[nnr_x,K,2] - 2Z[J,K,2] + Z[nnl_x,K,2])/dx^2 + (Z[nnr_y,K,2] - 2Z[J,K,2] + Z[nnl_y,K,2])/dy^2
                - ( m_ρ^2 + α*abs2(ϕ[J,2]) + β*ψ[J,2]^2 ) * Z[J,K,2] )
end

function updateForNextStep!(ϕ,ψ,Z,dϕdt,dψdt,dZdt)
    @views begin
        ϕ[:,1] .= ϕ[:,2]
        ψ[:,1] .= ψ[:,2]
        Z[:,:,1] .= Z[:,:,2]
        dϕdt[:,1] .= dϕdt[:,2]
        dψdt[:,1] .= dψdt[:,2]
        dZdt[:,:,1] .= dZdt[:,:,2]
    end
end

# function recordSnap(ϕ,ψ,totalE,ZED_t,
#                     ϕdataIO,ψdataIO,energyIO,ZedIO)

#     # writedlm(ϕdataIO, ϕ)
#     # writedlm(ψdataIO, ψ)
#     # writedlm(energyIO,totalE)
#     # writedlm(ZedIO,ZED_t)

#     # write(ϕJLD,"snapshot_$count",ϕ[:,:,0])#!
#     # write(ψJLD,"snapshot_$count",ψ[:,:,0])#!
# end

end #module





# #Crank-Nicolson Method
    # export time_evolve!
    # function time_evolve!(ϕ,ψ,Z,dϕdt,dψdt,dZdt,meanSqrRenorm,zPE)

    #     #Flatten ϕ, ψ, dϕdt, dψdt #!no offset for these#!
    #     ϕ_s = im*zeros(N^2,2)
    #     ψ_s = zeros(N^2,2)
    #     dϕdt_s = im*zeros(N^2,2)
    #     dψdt_s = zeros(N^2,2)
    #     ϕ_s[:,1] =  @views flattenDimension(ϕ[:,:,0])   #!for these .= uses more allocations!
    #     ψ_s[:,1] =  @views flattenDimension(ψ[:,:,0])  
    #     dϕdt_s[:,1] =  @views flattenDimension(dϕdt[:,:,0])
    #     dψdt_s[:,1] =  @views flattenDimension(dψdt[:,:,0])



    #     #Snapshotting interval
    #     if round(nt/nsnaps,RoundDown) == 0
    #         snapInterval = 1
    #         println("---Caution: Number of time steps is smaller than snaps!---")
    #         println("Instead of ",nsnaps," snapshots, ",nt," snaps will be taken.","---")
    #     else
    #         snapInterval = round(Int,nt/nsnaps,RoundDown)
    #         if nt/nsnaps != snapInterval
    #             actualSnaps=round(Int,nt/snapInterval,RoundDown)
    #             println("---Caution: Number of snapshots are not multiple of time steps.")
    #             println("instead of ",nsnaps," snapshots,",actualSnaps," snaps will be taken.","---")
    #         end
    #     end

    #     ϕdataIO = open("data/phi.dat","w") 
    #     ψdataIO = open("data/psi.dat","w")
    #     energyIO = open("data/energy.dat","w")
    #     ZedIO = open("data/energies/ZED.dat","w")
    #     meanRhoSqrIO = open("data/meanRhoSqr.dat","w")

    #     for t=1:nt
    #         #Move a time step               #t
    #         eulerTo_Intermediate!(ϕ_s,ψ_s,Z,dϕdt_s,dψdt_s,dZdt,meanSqrRenorm)   #tilde1
    #         averageHalfStep!(ϕ_s,ψ_s,Z,dϕdt_s,dψdt_s,dZdt)        #bar1
    #         leap_forward!(ϕ_s,ψ_s,Z,dϕdt_s,dψdt_s,dZdt,meanSqrRenorm)          #tilde2
    #         averageHalfStep!(ϕ_s,ψ_s,Z,dϕdt_s,dψdt_s,dZdt)        #bar2
    #         leap_forward!(ϕ_s,ψ_s,Z,dϕdt_s,dψdt_s,dZdt,meanSqrRenorm)          #t+1
    #         #Shift next time values to present time for the next step
    #         updateForNextStep!(ϕ_s,ψ_s,Z,dϕdt_s,dψdt_s,dZdt)

    #         # Take a snap
    #         if mod(t,snapInterval) == 0
    #             ϕ[:,:,0]    .= ravelDimension(ϕ_s[:,1])    
    #             ψ[:,:,0]    .= ravelDimension(ψ_s[:,1])
    #             dϕdt[:,:,0] .= ravelDimension(dϕdt_s[:,1])
    #             dψdt[:,:,0] .= ravelDimension(dψdt_s[:,1])
    #             writedlm(ϕdataIO,ϕ[:,:,0])
    #             writedlm(ψdataIO,ψ[:,:,0])

    #             totalE,ZED = energy(ϕ,ψ,Z,dϕdt,dψdt,dZdt,meanSqrRenorm,zPE)
    #             writedlm(energyIO,totalE)
    #             ZED_two = ravelDimension(ZED)
    #             writedlm(ZedIO,ZED_two)


    #             #2-point function
    #             meanSqr_Rho = zeros(N^2)
    #             for J=1:N^2
    #                 meanSqr_Rho[J] = sum(abs2, Z[J,:,1]) #- meanSqrRenorm
    #             end
    #             meanSqr_Rho_two =ravelDimension(meanSqr_Rho)
    #             writedlm(meanRhoSqrIO,meanSqr_Rho_two)

    #         end

    #     end

    # end

    # function eulerTo_Intermediate!(ϕ,ψ,Z,dϕdt,dψdt,dZdt,meanSqrRenorm)

    #     ϕ_flux, ψ_flux = fluxes_ϕ_ψ(ϕ[:,1],ψ[:,1],Z[:,:,1],meanSqrRenorm)
    #     Z_flux = flux_Z(ϕ[:,1],ψ[:,1],Z[:,:,1])

    #     ϕ[:,2]      .= ϕ[:,1] .+ dt.* ( dϕdt[:,1] )
    #     ψ[:,2]      .= ψ[:,1] .+ dt.* ( dψdt[:,1] )
    #     Z[:,:,2]    .= Z[:,:,1] .+ dt.* ( dZdt[:,:,1] )
    #     dϕdt[:,2]   .= dϕdt[:,1] .+ dt.* ( ϕ_flux )
    #     dψdt[:,2]   .= dψdt[:,1] .+ dt.* ( ψ_flux )
    #     dZdt[:,:,2] .= dZdt[:,:,1] .+ dt.*( Z_flux )

    #     # for J=1:N^2

    #     #     # ϕ_flux, ψ_flux = fluxes_ϕ_ψ(ϕ[:,1],ψ[:,1],Z[:,:,1],meanSqrRenorm,J,t)
    #     #     # Z_flux = flux_Z(ϕ[:,1],ψ[:,1],Z[:,:,1],J,t)
        
    #     #     dϕdt[J,2] = dϕdt[J,1] + dt*( ϕ_flux )
    #     #     dψdt[J,2] = dψdt[J,1] + dt*( ψ_flux )
    #     #     ϕ[J,2] = ϕ[J,1] + dt*( dϕdt[J,1] )
    #     #     ψ[J,2] = ψ[J,1] + dt*( dψdt[J,1] )
    #     #     for K=1:N^2
    #     #         dZdt[J,K,2] = dZdt[J,K,1] + dt*( Z_flux )
    #     #         Z[J,K,2] = Z[J,K,1] + dt*( dZdt[J,K,1] )
    #     #     end
    #     # end


    # end

    # function averageHalfStep!(ϕ,ψ,Z,dϕdt,dψdt,dZdt)
    #     ϕ[:,2]      .= 1/2 .* (ϕ[:,1] .+ ϕ[:,2])
    #     ψ[:,2]      .= 1/2 .* (ψ[:,1] .+ ψ[:,2])
    #     Z[:,:,2]    .= 1/2 .* (Z[:,:,1] .+ Z[:,:,2])
    #     dϕdt[:,2]   .= 1/2 .* (dϕdt[:,1] .+ dϕdt[:,2])
    #     dψdt[:,2]   .= 1/2 .* (dψdt[:,1] .+ dψdt[:,2])
    #     dZdt[:,:,2] .= 1/2 .* (dZdt[:,:,1] .+ dZdt[:,:,2])

    #     # for J=1:N^2
    #     #     dϕdt[J,2] = ( dϕdt[J,1] + dϕdt[J,2] )/2
    #     #     dψdt[J,2] = ( dψdt[J,1] + dψdt[J,2] )/2
    #     #     ϕ[J,2] = ( ϕ[J,1] + ϕ[J,2] )/2
    #     #     ψ[J,2] = ( ψ[J,1] + ψ[J,2] )/2
    #     #     for K=1:N^2
    #     #         dZdt[J,K,2] = ( dZdt[J,K,1] + dZdt[J,K,2] )/2
    #     #         Z[J,K,2] = ( Z[J,K,1] + Z[J,K,2] )/2
    #     #     end
    #     # end
    # end

    # function leap_forward!(ϕ,ψ,Z,dϕdt,dψdt,dZdt,meanSqrRenorm)

    #     ϕ_flux, ψ_flux = fluxes_ϕ_ψ(ϕ[:,2],ψ[:,2],Z[:,:,2],meanSqrRenorm)
    #     Z_flux = flux_Z(ϕ[:,2],ψ[:,2],Z[:,:,2])

    #     # ϕ[:,2]      .= ϕ[:,1] .+ dt.*( dϕdt[:,2] )
    #     # ψ[:,2]      .= ψ[:,1] .+ dt.*( dψdt[:,2] )
    #     # Z[:,:,2]    .= Z[:,:,1] .+ dt.*( dZdt[:,:,2] )  
    #     # dϕdt[:,2]   .= dϕdt[:,1] .+ dt.*( ϕ_flux )
    #     # dψdt[:,2]   .= dψdt[:,1] .+ dt.*( ψ_flux )
    #     # dZdt[:,:,2] .= dZdt[:,:,1] .+ dt.*( Z_flux )

    #     ϕ_temp    = ϕ[:,1] .+ dt.*( dϕdt[:,2] )
    #     ψ_temp    = ψ[:,1] .+ dt.*( dψdt[:,2] )
    #     Z_temp    = Z[:,:,1] .+ dt.*( dZdt[:,:,2] )  
    #     dϕdt_temp = dϕdt[:,1] .+ dt.*( ϕ_flux )
    #     dψdt_temp = dψdt[:,1] .+ dt.*( ψ_flux )
    #     dZdt_temp = dZdt[:,:,1] .+ dt.*( Z_flux )

    #     ϕ[:,2]      .= ϕ_temp
    #     ψ[:,2]      .= ψ_temp
    #     Z[:,:,2]    .= Z_temp
    #     dϕdt[:,2]   .= dϕdt_temp
    #     dψdt[:,2]   .= dψdt_temp
    #     dZdt[:,:,2] .= dZdt_temp

    # end

    # function fluxes_ϕ_ψ(ϕ,ψ,Z,meanSqrRenorm)#,ψ,Z,meanSqrRenorm,J,nnl_x,nnr_x,nnl_y,nnr_y)
    # # function fluxes_ϕ_ψ(ϕ,ψ,Z,meanSqrRenorm,J,t)

    #     ϕ_flux = zeros(N^2)
    #     ψ_flux = zeros(N^2)
    #     for J=1:N^2
    #         nnl_x, nnr_x, nnl_y, nnr_y = pbc1D(J)
    #         #2-point function
    #         meanSqr_Rho = sum(abs2, Z[J,:])

    #         ϕ_flux[J] = ( (ϕ[nnr_x] 
    #                         - 2ϕ[J] 
    #                         + ϕ[nnl_x])/dx^2 
    #                     + (ϕ[nnr_y] - 2ϕ[J] + ϕ[nnl_y])/dy^2 
    #                         + ( m_ϕ^2 -α/2 *(meanSqr_Rho - meanSqrRenorm)/(dx*dy) -λ* abs2(ϕ[J]) ) * ϕ[J] )
    #         ψ_flux[J] = ( (ψ[nnr_x] - 2ψ[J] + ψ[nnl_x])/dx^2 + (ψ[nnr_y] - 2ψ[J] + ψ[nnl_y])/dy^2
    #                         - ( m_ψ^2 +β *(meanSqr_Rho - meanSqrRenorm)/(dx*dy) ) * ψ[J] )
    #     end

    #     # ϕ_flux = zeros(N^2)
    #     # ψ_flux = zeros(N^2)
    #     # nnl_x, nnr_x, nnl_y, nnr_y = pbc1D(J)
    #     # #2-point function
    #     # meanSqr_Rho = 0.0
    #     # for K=1:N^2
    #     #     meanSqr_Rho = meanSqr_Rho + abs2(Z[J,K,t]) 
    #     # end

    #     # ϕ_flux[J] = ( (ϕ[nnr_x] 
    #     #                 - 2ϕ[J] 
    #     #                 + ϕ[nnl_x])/dx^2 
    #     #                 + (ϕ[nnr_y] - 2ϕ[J] + ϕ[nnl_y])/dy^2 
    #     #                 + ( m_ϕ^2 -α/2 *(meanSqr_Rho - meanSqrRenorm)/(dx*dy) -λ* abs2(ϕ[J]) ) * ϕ[J] )
    #     # ψ_flux[J] = ( (ψ[nnr_x] - 2ψ[J] + ψ[nnl_x])/dx^2 + (ψ[nnr_y] - 2ψ[J] + ψ[nnl_y])/dy^2
    #     #                     - ( m_ψ^2 +β *(meanSqr_Rho - meanSqrRenorm)/(dx*dy) ) * ψ[J] )
    #     # Z_flux = zeros(N^2,N^2)
    #     # for K=1:N^2
    #     #     Z_flux[J,K] = ( (Z[nnr_x,K,t] - 2Z[J,K,t] + Z[nnl_x,K,t])/dx^2 + (Z[nnr_y,K,t] - 2Z[J,K,t] + Z[nnl_y,K,t])/dy^2
    #     #                     - ( m_ρ^2 + α*abs2(ϕ[J]) + β*ψ[J]^2 ) * Z[J,K,t] )
    #     # end


    # return ϕ_flux, ψ_flux#, Z_flux
    # end

    # function flux_Z(ϕ,ψ,Z)
    #     Z_flux = im*zeros(N^2,N^2)
    #     for J=1:N^2
    #         nnl_x, nnr_x, nnl_y, nnr_y = pbc1D(J)
    #         for K=1:N^2
    #             Z_flux[J,K] = ( (Z[nnr_x,K] 
    #                             - 2Z[J,K] 
    #                             + Z[nnl_x,K])/dx^2 +
    #                              (Z[nnr_y,K] - 2Z[J,K] + Z[nnl_y,K])/dy^2
    #                     - ( m_ρ^2 + α*abs2(ϕ[J]) + β*ψ[J]^2 ) * Z[J,K] )
    #         end
    #     end

    # return Z_flux
    # end

    # function updateForNextStep!(ϕ,ψ,Z,dϕdt,dψdt,dZdt)
    #     ϕ[:,1] .= ϕ[:,2]
    #     ψ[:,1] .= ψ[:,2]
    #     Z[:,:,1] .= Z[:,:,2]
    #     dϕdt[:,1] .= dϕdt[:,2]
    #     dψdt[:,1] .= dψdt[:,2]
    #     dZdt[:,:,1] .= dZdt[:,:,2]
    # end 

    # end





# # 4-index Notation (Leap-frog)
# export time_evolve!
# function time_evolve!(ϕ,ψ,Z,dϕdt,dψdt,dZdt,meanSqrRenorm,zPE)

#     Z_t = Array{ComplexF64,5}(undef, Nx,Nx,Ny,Ny,2) #This way seems to be faster and memory friendly for complex arrays. 
#     dZdt_t = Array{ComplexF64,5}(undef, Nx,Nx,Ny,Ny,2)
#     Z_t = OffsetArray(Z_t,lx:rx,lx:rx,ly:ry,ly:ry,0:1)     
#     dZdt_t = OffsetArray(dZdt_t,lx:rx,lx:rx,ly:ry,ly:ry,0:1)
#     #Map Z to 4 index
#     Z_t[:,:,:,:,0] .= @views mapZTo4Index(Z[:,:,1])
#     dZdt_t[:,:,:,:,0] .= @views mapZTo4Index(dZdt[:,:,1])

#     #Snapshotting interval
#     if round(nt/nsnaps,RoundDown) == 0
#         snapInterval = 1
#         println("---Caution: Number of time steps is smaller than snaps!---")
#         println("Instead of ",nsnaps," snapshots, ",nt," snaps will be taken.","---")
#     else
#         snapInterval = round(Int,nt/nsnaps,RoundDown)
#         if nt/nsnaps != snapInterval
#             actualSnaps=round(Int,nt/snapInterval,RoundDown)
#             println("---Caution: Number of snapshots are not multiple of time steps.")
#             println("Instead of ",nsnaps," snapshots, ",actualSnaps," snaps will be taken.","---")
#         end
#     end

#     ϕdataIO = open("data/phi.dat","w") 
#     ψdataIO = open("data/psi.dat","w")
#     energyIO = open("data/energy.dat","w")
#     ZedIO = open("data/energies/ZED.dat","w")


#     for t=1:nt

#         #Move a time step
#         half_step!(ϕ,ψ,Z_t,dϕdt,dψdt,dZdt_t,0)
#         leap_forward!(ϕ,ψ,Z_t,dϕdt,dψdt,dZdt_t,meanSqrRenorm)
#         half_step!(ϕ,ψ,Z_t,dϕdt,dψdt,dZdt_t,1)
#         #Shift next time values to present time for the next step
#         updateForNextStep(ϕ,ψ,Z_t,dϕdt,dψdt,dZdt_t)

#         # #!map z back to 2 index for printing easy matrix notation
#         # Z[:,:,1] .= @views mapZTo2Index(Z_t[:,:,:,:,0])
#         # dZdt[:,:,1] .= @views mapZTo2Index(dZdt_t[:,:,:,:,0])

#         # Take a snap
#         if mod(t,snapInterval) == 0

#             # totalE,ZED = energy(ϕ,ψ,Z,dϕdt,dψdt,dZdt,meanSqrRenorm,zPE)
#             totalE,ZED = energy(ϕ,ψ,Z_t,dϕdt,dψdt,dZdt_t,meanSqrRenorm,zPE)
#             ZED_t = ravelDimension(ZED) #!unnecessary temp array.

#             writedlm(ϕdataIO, ϕ[:,:,0])
#             writedlm(ψdataIO, ψ[:,:,0])
#             writedlm(energyIO,totalE)
#             writedlm(ZedIO,ZED_t)
        

#             #!
#             # @views recordSnap(ϕ[:,:,0],ψ[:,:,0],totalE,ZED_t,
#             #                     ϕdataIO,ψdataIO,energyIO,ZedIO)
    
#             #!
#             #check constraints
#             # @views constraints_checker(Z[:,:,1],dZdt[:,:,1])
#             # @views conserved_checker(Z[:,:,1],dZdt[:,:,1])
#         end
#     end
# end


# function half_step!(ϕ,ψ,Z,dϕdt,dψdt,dZdt,t)
#     dt_half =dt/2   #!not sure yet if i wanna keep them. harder to read.
#     for j=lx:rx
#         for k=ly:ry
#             ϕ[j,k,1] = ϕ[j,k,t] + dt_half*( dϕdt[j,k,t] ) 
#             ψ[j,k,1] = ψ[j,k,t] + dt_half*( dψdt[j,k,t] )
#             for l=lx:rx
#                 for m=ly:ry
#                     Z[j,l,k,m,1] = Z[j,l,k,m,t] + dt_half*( dZdt[j,l,k,m,t] )
#                 end
#             end
#         end
#     end
# end


# function leap_forward!(ϕ,ψ,Z,dϕdt,dψdt,dZdt,meanSqrRenorm)

#     for j=lx:rx
#         for k=ly:ry
#             #PBC
#             nnl_x, nnr_x, nnl_y, nnr_y = pbc2D(j,k)
#             #Calculate fluxes for ϕ and ψ
#             @inline ϕ_flux , ψ_flux =  @views fluxes_ϕ_ψ(ϕ[:,:,1],ψ[:,:,1],Z[:,:,:,:,1],meanSqrRenorm,j,k,nnl_x,nnr_x,nnl_y,nnr_y)
#             dϕdt[j,k,1] = dϕdt[j,k,0] + dt*( ϕ_flux )
#             dψdt[j,k,1] = dψdt[j,k,0] + dt*( ψ_flux )
#             for l=lx:rx
#                 for m=ly:ry
#                     # @views Z_flux = flux_Z(ϕ[j,k,1],ψ[j,k,1],Z[:,l,:,m,1],j,k,nnl_x,nnr_x,nnl_y,nnr_y)
#                     @inline Z_flux = flux_Z(ϕ,ψ,Z,j,l,k,m,nnl_x,nnr_x,nnl_y,nnr_y)
#                     dZdt[j,l,k,m,1] = dZdt[j,l,k,m,0] + dt*( Z_flux )
#                 end
#             end
#         end
#     end

# end

# @inline function fluxes_ϕ_ψ(ϕ,ψ,Z,meanSqrRenorm,j,k,nnl_x,nnr_x,nnl_y,nnr_y)
#     #2-point function
#     meanSqr_Rho = @views sum(abs2, Z[j,:,k,:]) #!i need to check if this is the same as 2-index notation.

#     ϕ_flux = ( (ϕ[nnr_x,k] - 2ϕ[j,k] + ϕ[nnl_x,k])/dx^2 + (ϕ[j,nnr_y] - 2ϕ[j,k] + ϕ[j,nnl_y])/dy^2 
#             + ( m_ϕ^2 - α/2 *(meanSqr_Rho - meanSqrRenorm)/(dx*dy) -λ* abs2(ϕ[j,k]) ) * ϕ[j,k] )
#     ψ_flux = ( (ψ[nnr_x,k] - 2ψ[j,k] + ψ[nnl_x,k])/dx^2 + (ψ[j,nnr_y] - 2ψ[j,k] + ψ[j,nnl_y])/dy^2 
#             - ( m_ψ^2 + β *(meanSqr_Rho - meanSqrRenorm)/(dx*dy) ) * ψ[j,k] )

# return ϕ_flux ,ψ_flux
# end

# @inline function flux_Z(ϕ,ψ,Z,j,l,k,m,nnl_x,nnr_x,nnl_y,nnr_y)
#     Z_flux = ( (Z[nnr_x,l,k,m,1] - 2Z[j,l,k,m,1] + Z[nnl_x,l,k,m,1])/dx^2 
#              + (Z[j,l,nnr_y,m,1] - 2Z[j,l,k,m,1] + Z[j,l,nnl_y,m,1])/dy^2
#              - ( m_ρ^2 + α*abs2(ϕ[j,k,1]) + β*ψ[j,k,1]^2 ) * Z[j,l,k,m,1] )
# end

# function updateForNextStep(ϕ,ψ,Z,dϕdt,dψdt,dZdt)
#     @views begin
#         ϕ[:,:,0] .= ϕ[:,:,1]
#         ψ[:,:,0] .= ψ[:,:,1]
#         Z[:,:,:,:,0] .= Z[:,:,:,:,1]
#         dϕdt[:,:,0] .= dϕdt[:,:,1]
#         dψdt[:,:,0] .= dψdt[:,:,1]
#         dZdt[:,:,:,:,0] .= dZdt[:,:,:,:,1]
#     end
# end

# end #module
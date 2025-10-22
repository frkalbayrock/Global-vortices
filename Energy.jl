module Energy

using Distributed
include("Parameters.jl")
using .Parameters
include("Auxiliary.jl")
using .Auxiliary_Routines



export energy!
function energy!(ZED_gl,meanSqrRenorm,zPE)

    #--Update the paddings before calculating energy
    @everywhere workers() update_Paddings!(ϕ,ψ,Z)

    #--Calculate energy on each chunk
    @everywhere workers() totalE_loc = energy_calculation!(ϕ,ψ,Z,dϕdt,dψdt,dZdt,ZED,($meanSqrRenorm),($zPE))


    #--Combine local energies for total energy 
    totalE = 0.
    for i in workers()
        totalE += (@fetchfrom i Main.totalE_loc)::Float64
    end
    

    #---Collect ZED to global
    for i in workers()
        lx_p, rx_p, ly_p, ry_p = chunker(i)
        ZED_gl[lx_p:rx_p, ly_p:ry_p] .= (@fetchfrom i Main.ZED[:,:])::Array{Float64, 2}
    end

return totalE
end




@everywhere workers() function energy_calculation!(ϕ,ψ,Z,dϕdt,dψdt,dZdt,ZED,meanSqrRenorm,zPE)

    #---Using Two-Index Notation---#

        #Classical fields densities
        kEϕψ = 0.
        gEϕψ = 0.
        pEϕψ = 0.
        classicalE = 0.
        #Quantum field density
        quantumEZ = 0.
        #Total energy density (classical+quantum)
        sumE = 0.
        #Total Energies
        totalE = 0.


        #Integrate over all space
        for j=padd+1:padd+Nx_loc
            for k=padd+1:padd+Ny_loc
                # #Get neighbours with pbc
                # nnl_x, nnr_x, nnl_y, nnr_y = pbc2D(j,k)

                #Kinetic, gradient, potential energy densities
                kEϕψ = abs2(dϕdt[j-padd,k-padd,1]) + (dψdt[j-padd,k-padd,1])^2/2
                gEϕψ = ( ( ( (ψ[j+1,k] - ψ[j,k])/dx )^2  
                         + ( (ψ[j,k+1] - ψ[j,k])/dy )^2 )/2
                        + ( abs2( (ϕ[j+1,k] - ϕ[j,k])/dx ) 
                          + abs2( (ϕ[j,k+1] - ϕ[j,k])/dy ) ) )
                pEϕψ = -m_ϕ^2*abs2(ϕ[j,k]) + m_ψ^2*ψ[j,k]^2/2 + λ*abs2(ϕ[j,k])^2/2 + λ*η^4/2
                #Total classical energy density
                classicalE = kEϕψ + gEϕψ + pEϕψ

                #Energy densities of Z
                kEZ = 0.
                gEZ = gEZfwd = gEZbkd = 0.
                pEZ = 0.
                intEZ = 0.

                for l=1:Nx
                    for m=1:Ny
                        kEZ = kEZ + ( abs2(dZdt[j-padd,l,k-padd,m,1]) )/2
                        gEZfwd = ( abs2((Z[j+1,l,k,m] - Z[j,l,k,m])/dx)  
                                 + abs2((Z[j,l,k+1,m] - Z[j,l,k,m])/dy) )/2
                        gEZbkd = ( abs2((Z[j,l,k,m] - Z[j-1,l,k,m])/dx)  
                                 + abs2((Z[j,l,k,m] - Z[j,l,k-1,m])/dy) )/2
                        gEZ = gEZ + (gEZbkd + gEZfwd)/2
                        pEZ = pEZ + m_ρ^2* abs2(Z[j,l,k,m])/2
                        intEZ = intEZ + ( α*abs2(ϕ[j,k]) + β*ψ[j,k]^2 ) * ( abs2(Z[j,l,k,m]) )/2
                    end
                end
            
                #Renormalization
                intEZ = intEZ - ( α*abs2(ϕ[j,k]) + β*ψ[j,k]^2 )*meanSqrRenorm/2
                #Total energy density of Z including the interactions
                quantumEZ = (kEZ + gEZ + pEZ + intEZ)/(dx*dy)

                #Zero-Point Energy
                quantumEZ = quantumEZ - zPE

                #Total energy density
                sumE = classicalE + quantumEZ

                #Total energy
                totalE = totalE + dx*dy*sumE

                #Energy Density of Z field
                ZED[j-padd,k-padd] = quantumEZ

            end
        end

return totalE
end




end #module
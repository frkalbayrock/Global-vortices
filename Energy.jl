module Energy

using Distributed
include("Parameters.jl")
using .Parameters
include("Auxiliary.jl")
using .Auxiliary_Routines
using OffsetArrays
using DistributedArrays
#!
@everywhere using Profile, PProf, InteractiveUtils

export energy
function energy(ϕ,ψ,Z,dϕdt,dψdt,dZdt,ZED,meanSqrRenorm,zPE)


    # #!
    # @sync @distributed for _ in workers()
    #     lan = size(localpart(ψ))
    #     anan = zeros(lan[1],lan[2]).+myid()
    #     @views localpart(ψ)[:,:,1] .= anan
        
    #     anan = im*zeros(lan[1],lan[2]).+myid()
    #     @views localpart(ϕ)[:,:,1] .= anan

    # end
    # #!

    #--Update the paddings before calculating energy
    update_Paddings!(ϕ,ψ,Z,1)


    #--Calculate the total energy
    totalE = @sync @distributed (+) for _ in workers()
        totalE_loc = energy_calculation!(localpart(ϕ),localpart(ψ),localpart(Z),
                                         localpart(dϕdt),localpart(dψdt),localpart(dZdt),
                                         localpart(ZED),meanSqrRenorm,zPE)
        totalE_loc
    end


return totalE
end


function energy_calculation!(ϕ,ψ,Z,dϕdt,dψdt,dZdt,ZED,meanSqrRenorm,zPE)


    #---Using Two-Index Notation---#
        # Z_t , dZdt_t = mapZTo4Index(Z[:,:,1],dZdt[:,:,1]) #!This needs an update since the definition of the function changed.
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
                gEϕψ = ( ( ( (ψ[j+1,k,1] - ψ[j,k,1])/dx )^2  
                         + ( (ψ[j,k+1,1] - ψ[j,k,1])/dy )^2 )/2
                        + ( abs2( (ϕ[j+1,k,1] - ϕ[j,k,1])/dx ) 
                          + abs2( (ϕ[j,k+1,1] - ϕ[j,k,1])/dy ) ) )
                pEϕψ = -m_ϕ^2*abs2(ϕ[j,k,1]) + m_ψ^2*ψ[j,k,1]^2/2 + λ*abs2(ϕ[j,k,1])^2/2 + λ*η^4/2
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
                        gEZfwd = ( abs2((Z[j+1,l,k,m,1] - Z[j,l,k,m,1])/dx)  
                                 + abs2((Z[j,l,k+1,m,1] - Z[j,l,k,m,1])/dy) )/2
                        gEZbkd = ( abs2((Z[j,l,k,m,1] - Z[j-1,l,k,m,1])/dx)  
                                 + abs2((Z[j,l,k,m,1] - Z[j,l,k-1,m,1])/dy) )/2
                        gEZ = gEZ + (gEZbkd + gEZfwd)/2
                        pEZ = pEZ + m_ρ^2* abs2(Z[j,l,k,m,1])/2
                        intEZ = intEZ + ( α*abs2(ϕ[j,k,1]) + β*ψ[j,k,1]^2 ) * ( abs2(Z[j,l,k,m,1]) )/2
                    end
                end
            
                #Renormalization
                intEZ = intEZ - ( α*abs2(ϕ[j,k,1]) + β*ψ[j,k,1]^2 )*meanSqrRenorm/2
                #Total energy density of Z including the interactions
                quantumEZ = (kEZ + gEZ + pEZ + intEZ)/(dx*dy)

                #Zero-Point Energy
                quantumEZ = quantumEZ - zPE

                #Total energy density
                sumE = classicalE + quantumEZ

                #Total energy
                totalE = totalE + dx*dy*sumE

                #Energy density of Z update
                ZED[j-padd,k-padd] = quantumEZ

            end
        end

        
return totalE
end


# #---Using Single-Index Notation---#   
        #     #2-index to 1-index mapping
        #     ϕ_s = @views flattenDimension(ϕ[:,:,0])
        #     ψ_s = @views flattenDimension(ψ[:,:,0])
        #     dϕdt_s = @views flattenDimension(dϕdt[:,:,0])
        #     dψdt_s = @views flattenDimension(dψdt[:,:,0])

        #     #Classical fields densities
        #     kEϕψ = 0.
        #     gEϕψ = 0.
        #     pEϕψ = 0.
        #     classicalE = 0.
        #     #Quantum field density
        #     quantumEZ = 0.
        #     #Total energy density (classical+quantum)
        #     sumE = 0.
        #     #Total Energies
        #     totalE = 0.

        #     #!
        #     ZED = zeros(1:N^2)
        #     # Edensity = zeros(N^2)

        #     #Integrate over all space
        #     for J=1:N^2
        #         #Get the neighbours with pbc
        #         nnl_x, nnr_x, nnl_y, nnr_y = pbc1D(J)

        #         #Kinetic(kE), gradient(gE) and potential(pE) energy densities
        #         kEϕψ = abs2(dϕdt_s[J]) + (dψdt_s[J])^2/2
        #         gEϕψfwd = ( ( ( (ψ_s[nnr_x] - ψ_s[J])/dx )^2  + ( (ψ_s[nnr_y] - ψ_s[J])/dy )^2 )/2
        #                     + ( abs2( (ϕ_s[nnr_x] - ϕ_s[J])/dx ) + abs2( (ϕ_s[nnr_y] - ϕ_s[J])/dy ) ) )
        #         gEϕψbkd = ( ( ( (ψ_s[J] - ψ_s[nnl_x])/dx )^2  + ( (ψ_s[J] - ψ_s[nnl_y])/dy )^2 )/2
        #                     + ( abs2( (ϕ_s[J] - ϕ_s[nnl_x])/dx ) + abs2( (ϕ_s[J] - ϕ_s[nnl_y])/dy ) ) )
        #         gEϕψ = (gEϕψfwd + gEϕψbkd)/2
        #         pEϕψ = -m_ϕ^2*abs2(ϕ_s[J]) + m_ψ^2*ψ_s[J]^2/2 + λ*abs2(ϕ_s[J])^2/2 #+ λ*η^4/2
        #         #Total classical energy density
        #         classicalE = kEϕψ + gEϕψ + pEϕψ

        #         #Energy densities of Z
        #         kEZ = 0.
        #         gEZ = gEZfwd = gEZbkd = 0.
        #         pEZ = 0.
        #         intEZ = 0.
        #         for K=1:N^2
        #             kEZ = kEZ + ( abs2(dZdt[J,K,1]) )/2
        #             gEZfwd = ( abs2((Z[nnr_x,K,1] - Z[J,K,1])/dx)  
        #                      + abs2((Z[nnr_y,K,1] - Z[J,K,1])/dy) )/2
        #             gEZbkd = ( abs2((Z[J,K,1] - Z[nnl_x,K,1])/dx)  
        #                      + abs2((Z[J,K,1] - Z[nnl_y,K,1])/dy) )/2
        #             gEZ = gEZ + (gEZbkd + gEZfwd)/2
        #             pEZ = pEZ + m_ρ^2* abs2(Z[J,K,1])/2
        #             intEZ = intEZ + ( α*abs2(ϕ_s[J]) + β*ψ_s[J]^2 ) * ( abs2(Z[J,K,1]) )/2
        #         end

        #         #-Renormalization 
        #         intEZ = intEZ - (α*abs2(ϕ_s[J]) + β*ψ_s[J]^2)*meanSqrRenorm/2
        #         #Total energy density of Z including the interactions
        #         quantumEZ = (kEZ + gEZ + pEZ + intEZ)/(dx*dy)

        #         #Zero-Point Energy
        #         quantumEZ = quantumEZ - zPE

        #         #-Total energy density
        #         sumE =  classicalE + quantumEZ

        #         #-Total energy integral
        #         totalE = totalE + dx*dy*sumE

        #         #!Change later
        #         ZED[J] = quantumEZ
        #         # Edensity[J] = classicalE
        #     end

        #     #!Change later
        #     # return ZED
        #     return totalE,ZED#, Edensity
#

end #module
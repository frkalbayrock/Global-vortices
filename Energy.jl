module Energy

include("Parameters.jl")
using .Parameters
using MPI
include("MPIAux.jl")
using .MPIAux
using DelimitedFiles

const comm = MPI.COMM_WORLD
const myrank = MPI.Comm_rank(comm)
const nprocs = MPI.Comm_size(comm)


export energy 
function energy(ϕ,ψ,Z,dϕdt,dψdt,dZdt,ZED,meanSqrRenorm,zPE)

    #--Update the paddings before calculating energy
    update_Paddings!(ϕ,ψ,Z)

    #--Calculate energy on each chunk
    totalE_loc, ZED_loc = energy_calculation(ϕ,ψ,Z,dϕdt,dψdt,dZdt,meanSqrRenorm,zPE)

    #--Combine local energies for total energy
    totalE = MPI.Reduce(totalE_loc,MPI.SUM,comm;root=0)
    if myrank!=0
        totalE=0
    end
    totalE = MPI.Bcast(totalE, 0, comm)

    #--Collect ZED to global
    recvbuff = @views MPI.gather(ZED_loc[1:Nx_loc,1:Ny_loc], comm; root=0)
    if myrank==0
        for rank=0:nprocs-1
            lx_loc , rx_loc, ly_loc, ry_loc = chunker(rank)
            ZED_loc = recvbuff[rank+1]
            Main.ZED_gl[lx_loc:rx_loc,ly_loc:ry_loc] .= @views ZED_loc[:,:]
        end
    end

return totalE
end


function energy_calculation(ϕ,ψ,Z,dϕdt,dψdt,dZdt,meanSqrRenorm,zPE)

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

    #!
    ZED_loc = zeros(Nx_loc,Ny_loc)

    #Integrate over all space
    for j=padd+1:padd+Nx_loc
        for k=padd+1:padd+Ny_loc

            #Kinetic, gradient, potential energy densities
            kEϕψ = abs2(dϕdt[j,k,1]) + (dψdt[j,k,1])^2/2
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
                    kEZ = kEZ + ( abs2(dZdt[j,l,k,m,1]) )/2
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

            #!Change later
            ZED_loc[j-padd,k-padd] = quantumEZ
        end
    end

    #!Change later
    return totalE,ZED_loc


# #---Using Single-Index Notation---#   
    # #2-index to 1-index mapping
    # ϕ_s = @views flattenDimension(ϕ[:,:,0])
    # ψ_s = @views flattenDimension(ψ[:,:,0])
    # dϕdt_s = @views flattenDimension(dϕdt[:,:,0])
    # dψdt_s = @views flattenDimension(dψdt[:,:,0])

    # #Classical fields densities
    # kEϕψ = 0.
    # gEϕψ = 0.
    # pEϕψ = 0.
    # classicalE = 0.
    # #Quantum field density
    # quantumEZ = 0.
    # #Total energy density (classical+quantum)
    # sumE = 0.
    # #Total Energies
    # totalE = 0.

    # #!
    # ZED = zeros(1:N^2)
    # # Edensity = zeros(N^2)

    # #Integrate over all space
    # for J=1:N^2
    #     #Get the neighbours with pbc
    #     nnl_x, nnr_x, nnl_y, nnr_y = pbc1D(J)

    #     #Kinetic(kE), gradient(gE) and potential(pE) energy densities
    #     kEϕψ = abs2(dϕdt_s[J]) + (dψdt_s[J])^2/2
    #     gEϕψfwd = ( ( ( (ψ_s[nnr_x] - ψ_s[J])/dx )^2  + ( (ψ_s[nnr_y] - ψ_s[J])/dy )^2 )/2
    #                 + ( abs2( (ϕ_s[nnr_x] - ϕ_s[J])/dx ) + abs2( (ϕ_s[nnr_y] - ϕ_s[J])/dy ) ) )

    #     gEϕψbkd = ( ( ( (ψ_s[J] - ψ_s[nnl_x])/dx )^2  + ( (ψ_s[J] - ψ_s[nnl_y])/dy )^2 )/2
    #                 + ( abs2( (ϕ_s[J] - ϕ_s[nnl_x])/dx ) + abs2( (ϕ_s[J] - ϕ_s[nnl_y])/dy ) ) )

    #     gEϕψ = (gEϕψfwd + gEϕψbkd)/2
    #     pEϕψ = -m_ϕ^2*abs2(ϕ_s[J]) + m_ψ^2*ψ_s[J]^2/2 + λ*abs2(ϕ_s[J])^2/2 #+ λ*η^4/2
    #     #Total classical energy density
    #     classicalE = kEϕψ + gEϕψ + pEϕψ

    #     #Energy densities of Z
    #     kEZ = 0.
    #     gEZ = gEZfwd = gEZbkd = 0.
    #     pEZ = 0.
    #     intEZ = 0.
    #     for K=1:N^2
    #         kEZ = kEZ + ( abs2(dZdt[J,K,1]) )/2
    #         gEZfwd = ( abs2((Z[nnr_x,K,1] - Z[J,K,1])/dx)  
    #                     + abs2((Z[nnr_y,K,1] - Z[J,K,1])/dy) )/2
    #         gEZbkd = ( abs2((Z[J,K,1] - Z[nnl_x,K,1])/dx)  
    #                     + abs2((Z[J,K,1] - Z[nnl_y,K,1])/dy) )/2
    #         gEZ = gEZ + (gEZbkd + gEZfwd)/2
    #         pEZ = pEZ + m_ρ^2* abs2(Z[J,K,1])/2
    #         intEZ = intEZ + ( α*abs2(ϕ_s[J]) + β*ψ_s[J]^2 ) * ( abs2(Z[J,K,1]) )/2
    #     end

    #     # #!
    #     # # Precompute reciprocal of dx and dy
    #     # inv_dx = 1 / dx
    #     # inv_dy = 1 / dy

    #     # # Initialize energy densities
    #     # kEZ = 0.0
    #     # gEZ = 0.0
    #     # pEZ = 0.0
    #     # intEZ = 0.0

    #     # # Loop over all elements
    #     # for K in 1:N^2
    #     #     zJK = Z[J, K, 1]
            
    #     #     kEZ += abs2(dZdt[J, K, 1]) / 2

    #     #     z_nnr_x = Z[nnr_x, K, 1]
    #     #     z_nnr_y = Z[nnr_y, K, 1]
    #     #     z_nnl_x = Z[nnl_x, K, 1]
    #     #     z_nnl_y = Z[nnl_y, K, 1]

    #     #     gEZfwd = (abs2((z_nnr_x - zJK) * inv_dx) + abs2((z_nnr_y - zJK) * inv_dy)) / 2
    #     #     gEZbkd = (abs2((zJK - z_nnl_x) * inv_dx) + abs2((zJK - z_nnl_y) * inv_dy)) / 2

    #     #     gEZ += (gEZbkd + gEZfwd) / 2

    #     #     pEZ += m_ρ^2 * abs2(zJK) / 2

    #     #     intEZ += (α * abs2(ϕ_s[J]) + β * ψ_s[J]^2) * abs2(zJK) / 2
    #     # end
    #     # #!

    #     #!
    #     # # Precompute reciprocal of dx and dy
    #     # inv_dx = 1 / dx
    #     # inv_dy = 1 / dy

    #     # # Initialize energy densities
    #     # kEZ = 0.0
    #     # gEZ = 0.0
    #     # pEZ = 0.0
    #     # intEZ = 0.0

    #     # # Preallocate arrays for intermediate calculations
    #     # gEZfwd = zeros(N^2)
    #     # gEZbkd = zeros(N^2)

    #     # # Loop over all elements
    #     # for K in 1:N^2
    #     #     zJK = Z[J, K, 1]

    #     #     kEZ += abs2(dZdt[J, K, 1]) / 2

    #     #     z_nnr_x = Z[nnr_x, K, 1]
    #     #     z_nnr_y = Z[nnr_y, K, 1]
    #     #     z_nnl_x = Z[nnl_x, K, 1]
    #     #     z_nnl_y = Z[nnl_y, K, 1]

    #     #     gEZfwd[K] = (abs2((z_nnr_x - zJK) * inv_dx) + abs2((z_nnr_y - zJK) * inv_dy)) / 2
    #     #     gEZbkd[K] = (abs2((zJK - z_nnl_x) * inv_dx) + abs2((zJK - z_nnl_y) * inv_dy)) / 2

    #     #     gEZ += (gEZbkd[K] + gEZfwd[K]) / 2

    #     #     pEZ += m_ρ^2 * abs2(zJK) / 2

    #     #     intEZ += (α * abs2(ϕ_s[J]) + β * ψ_s[J]^2) * abs2(zJK) / 2
    #     # end
    #     # #!



    #     #-Renormalization 
    #     intEZ = intEZ - (α*abs2(ϕ_s[J]) + β*ψ_s[J]^2)*meanSqrRenorm/2
    #     #Total energy density of Z including the interactions
    #     quantumEZ = (kEZ + gEZ + pEZ + intEZ)/(dx*dy)

    #     #Zero-Point Energy
    #     quantumEZ = quantumEZ - zPE

    #     #-Total energy density
    #     sumE =  classicalE + quantumEZ

    #     #-Total energy integral
    #     totalE = totalE + dx*dy*sumE

    #     #!Change later
    #     ZED[J] = quantumEZ
    #     # Edensity[J] = classicalE
    # end

    # #!Change later
    # # return ZED
    # return totalE,ZED#, Edensity
#

end

end #module
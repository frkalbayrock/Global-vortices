module Auxiliary_Routines 

include("Parameters.jl")
using .Parameters
include("IndexMap.jl")
using .IndexMap
using OffsetArrays

#This routine applies the periodic boundary conditions on 2-D lattice
#It takes the lattice point (j,k) and determines nnl (nearest neighbour left) and nnr (nearest neighbour right) on each direction.
#If the point (j,k) is at the boundary the nnl_x,y and nnr_x,y are set accordingly.
#Then it returns nnl_x,y and nnr_x,y.
export pbc2D
function pbc2D(j,k)

    #x-direction
    nnl_x = j-1
    nnr_x = j+1
    nnl_x < lx ?  nnl_x=rx : nothing
    nnr_x > rx ?  nnr_x=lx : nothing

    #y-direction
    nnl_y = k-1
    nnr_y = k+1
    nnl_y < ly ?  nnl_y=ry : nothing
    nnr_y > ry ?  nnr_y=ly : nothing

    return nnl_x, nnr_x, nnl_y, nnr_y 
end


#This routine find the nearest neighbours in x and y direction of a 2D lattice with pbc that is flattened. 
#It is assumed that the flattening is done using row-major order.
export pbc1D
function pbc1D(J)

    #-In x-direction:
    #   nearest neighbour on the right will be  +N away and
    #   nearest neighbour on the left will be -N away
    #   all under mod(N^2)
    if mod(J-N,N^2)!=0      #To make sure N^2 is mapped to N^2 not 0
        nnl_x = mod(J-N,N^2)    
    else
        nnl_x = N^2
    end

    if mod(J+N,N^2) !=0     #To make sure N^2 is mapped to N^2 not 0
        nnr_x = mod(J+N,N^2)
    else
        nnr_x = N^2
    end


    #-In y-direction:
    #   nearest neighbour on the right will be  +1 away and
    #   nearest neighbour on the left will be -1 away
    #   unless they are the last or first element of their block respectively
    #   For those we wrap around inside the block: 
    #   first to last (for the left neighbour) and last to first (for the right neighbour) is mapped.
    if mod(J,N) == 1
        nnl_y = J+N-1
    else
        nnl_y = J-1
    end

    if  mod(J,N) == 0
        nnr_y = J-N+1
    else
        nnr_y = J+1
    end

    return nnl_x, nnr_x, nnl_y, nnr_y
end

#
export energy
function energy(ϕ,ψ,Z,dϕdt,dψdt,dZdt,meanSqrRenorm,zPE)

    #---Using Single-Index Notation---#   
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
        #     pEϕψ = -m_ϕ^2*abs2(ϕ_s[J]) + m_ψ^2*ψ_s[J]^2/2 + λ*abs2(ϕ_s[J])^2/2 #+ λ*η^4/4
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
        #                  + abs2((Z[nnr_y,K,1] - Z[J,K,1])/dy) )/2
        #         gEZbkd = ( abs2((Z[J,K,1] - Z[nnl_x,K,1])/dx)  
        #                  + abs2((Z[J,K,1] - Z[nnl_y,K,1])/dy) )/2
        #         gEZ = gEZ + (gEZbkd + gEZfwd)/2
        #         pEZ = pEZ + m_ρ^2* abs2(Z[J,K,1])/2
        #         intEZ = intEZ + ( α*abs2(ϕ_s[J]) + β*ψ_s[J]^2 ) * ( abs2(Z[J,K,1]) )/2
        #     end

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
        ZED = zeros(Nx,Ny)
        ZED = OffsetArray(ZED,lx:rx,ly:ry)

        #Integrate over all space
        for j=lx:rx
            for k=ly:ry
                #Get neighbours with pbc
                nnl_x, nnr_x, nnl_y, nnr_y = pbc2D(j,k)

                #Kinetic, gradient, potential energy densities
                kEϕψ = abs2(dϕdt[j,k,0]) + (dψdt[j,k,0])^2/2
                gEϕψ = ( ( ( (ψ[nnr_x,k,0] - ψ[j,k,0])/dx )^2  + ( (ψ[j,nnr_y,0] - ψ[j,k,0])/dy )^2 )/2
                        + ( abs2( (ϕ[nnr_x,k,0] - ϕ[j,k,0])/dx ) + abs2( (ϕ[j,nnr_y,0] - ϕ[j,k,0])/dy ) ) )
                pEϕψ = -m_ϕ^2*abs2(ϕ[j,k,0]) + m_ψ^2*ψ[j,k,0]^2/2 + λ*abs2(ϕ[j,k,0])^2/2 #+ λ*eta^4/4
                #Total classical energy density
                classicalE = kEϕψ + gEϕψ + pEϕψ

                #Energy densities of Z
                kEZ = 0.
                gEZ = gEZfwd = gEZbkd = 0.
                pEZ = 0.
                intEZ = 0.

                for l=lx:rx
                    for m=ly:ry
                        kEZ = kEZ + ( abs2(dZdt[j,l,k,m,0]) )/2
                        gEZfwd = ( abs2((Z[nnr_x,l,k,m,0] - Z[j,l,k,m,0])/dx)  
                                    + abs2((Z[j,l,nnr_y,m,0] - Z[j,l,k,m,0])/dy) )/2
                        gEZbkd = ( abs2((Z[j,l,k,m,0] - Z[nnl_x,l,k,m,0])/dx)  
                                    + abs2((Z[j,l,k,m,0] - Z[j,l,nnl_y,m,0])/dy) )/2
                        gEZ = gEZ + (gEZbkd + gEZfwd)/2
                        pEZ = pEZ + m_ρ^2* abs2(Z[j,l,k,m,0])/2
                        intEZ = intEZ + ( α*abs2(ϕ[j,k,0]) + β*ψ[j,k,0]^2 ) * ( abs2(Z[j,l,k,m,0]) )/2
                    end
                end
            
                #Renormalization
                intEZ = intEZ - ( α*abs2(ϕ[j,k,0]) + β*ψ[j,k,0]^2 )*meanSqrRenorm/2
                #Total energy density of Z including the interactions
                quantumEZ = (kEZ + gEZ + pEZ + intEZ)/(dx*dy)

                #Zero-Point Energy
                quantumEZ = quantumEZ - zPE

                #Total energy density
                sumE = classicalE + quantumEZ

                #Total energy
                totalE = totalE + dx*dy*sumE

                #!Change later
                ZED[j,k] = quantumEZ
            end
        end

        #!Change later
        return totalE,ZED

end

end
module IC

include("Parameters.jl")
using .Parameters
using LinearAlgebra
include("IndexMap.jl")
using .IndexMap
using MPI
include("MPIAux.jl")
using .MPIAux

#!
using OffsetArrays
using DelimitedFiles

#!MPI PART
const comm = MPI.COMM_WORLD
const myrank = MPI.Comm_rank(comm)
const nprocs = MPI.Comm_size(comm)
const comm_cart = MPI.Cart_create(comm, nprocs_perdim; periodic=periods, reorder=false)
const coords_cart = MPI.Cart_coords(comm_cart, myrank)
const lx_p, rx_p, ly_p, ry_p = chunker(myrank) #p stands for "physical"


#initialConditions sets the initial conditions for the fields ϕ,ψ,Z and calculates the renormalization factor
export initialConditions!
#Set initial conditions for the fields ϕ, ψ, ρ
function initialConditions!(ϕ,ψ,Z,dϕdt,dψdt,dZdt)

    ic_ϕ_ψ!(ϕ,ψ,dϕdt,dψdt)

    #!For now let's just calculate the omega matrix on a single core. 
    recvbuff_ψ = @views MPI.gather(ψ[padd+1:padd+Nx_loc,padd+1:padd+Ny_loc], comm; root=0)
    recvbuff_ϕ = @views MPI.gather(ϕ[padd+1:padd+Nx_loc,padd+1:padd+Ny_loc], comm; root=0)
    #Reduce ϕ and ψ on master for calculation of Z
    if myrank==0
        for rank=0:nprocs-1
            lx_loc , rx_loc, ly_loc, ry_loc = chunker(rank)
            ϕ_loc = recvbuff_ϕ[rank+1]
            ψ_loc = recvbuff_ψ[rank+1]

            Main.ϕ_gl[lx_loc:rx_loc,ly_loc:ry_loc] .= ϕ_loc[:,:]
            Main.ψ_gl[lx_loc:rx_loc,ly_loc:ry_loc] .= ψ_loc[:,:]
        end

        #2-index to 1-index mapping
        ϕ_s = flattenDimension(Main.ϕ_gl)
        ψ_s = flattenDimension(Main.ψ_gl)
        #Get Sqrt(Omega) and its inverse matrices
        S_Ωzero, inv_S_Ωzero = omegaIC(ϕ_s,ψ_s)

        #!
        #Test 4-indexed S_Ωzero, inv_S_Ωzero
        S_Ωzero_f = mapZTo4Index(S_Ωzero)
        inv_S_Ωzero_f = mapZTo4Index(inv_S_Ωzero)

    end


    #Distribute Z for domain decomposition
    if myrank==0
        @views Z[lx_l:rx_l,:,ly_l:ry_l,:]    = -im/sqrt(2) .* inv_S_Ωzero_f[lx_p:rx_p,:,ly_p:ry_p,:]
        @views dZdt[lx_l:rx_l,:,ly_l:ry_l,:,1] =   1/sqrt(2) .* S_Ωzero_f[lx_p:rx_p,:,ly_p:ry_p,:]

        #---Send Z
        for idp=1:nprocs-1
             #Find the chunks location in the global
            lx_loc, rx_loc, ly_loc, ry_loc = chunker(idp)
            MPI.Send(inv_S_Ωzero_f[lx_loc:rx_loc,:,ly_loc:ry_loc,:].parent, comm; dest=idp, tag=1)
        end
        #---Send dZdt
        for idp=1:nprocs-1
            #Find the chunks location in the global
            lx_loc, rx_loc, ly_loc, ry_loc = chunker(idp)
            MPI.Send(S_Ωzero_f[lx_loc:rx_loc,:,ly_loc:ry_loc,:].parent, comm; dest=idp, tag=2)
        end
        
    else
        recvdata = zeros(Nx_loc,Nx,Ny_loc,Ny)
        MPI.Recv!(recvdata,comm; source=0,tag=1)
        Z[lx_l:rx_l,:,
          ly_l:ry_l,:] .= -im/sqrt(2) .* recvdata
        MPI.Recv!(recvdata,comm; source=0,tag=2)
        dZdt[lx_l:rx_l,:,
             ly_l:ry_l,:,1] .= 1/sqrt(2) .* recvdata
    end


end


#Sets the initial conditions for ϕ and ψ
function ic_ϕ_ψ!(ϕ,ψ,dϕdt,dψdt)
    width=2.0
    amp=10
    # vel=0.5
    vx=0.3
    vy=0.4
    r0=1.25
    γ=1/sqrt(1-(vx^2+vy^2))
    #ϕ and ψ i.c.
    for j in padd+1:Nx_loc+padd
        x = (lx_p+(j-padd)-1)*dx
        x1 = x-r0
        x2 = x+r0
        for k in padd+1:Ny_loc+padd
            y = (ly_p+(k-padd)-1)*dy
            y1 = y-r0
            y2 = y+r0
            ϕ[j,k] = η
            dϕdt[j,k,1] = 0
            ψ[j,k] = amp*(exp( -width/(vx^2 + vy^2)
                                * ( (x1 *(-vy) - y1 *(-vx))^2 + (x1* (-vx) + y1*(-vy))^2 * γ^2 ) )
                        + exp( -width/(vx^2 + vy^2) 
                                * ( (x2 *vy - y2 *vx)^2 + (x2* vx + y2 *vy)^2 * γ^2 )) )
            dψdt[j,k,1] = ( 2amp *width *γ^2 *(x1 *(-vx) + y1 *(-vy)) 
                            *exp(-width/(vx^2 + vy^2) * ( (x1 *(-vy) - y1 *(-vx))^2 + (x1* (-vx) + y1 *(-vy))^2 * γ^2 ))
                        + 2amp *width *γ^2 *(x2 *(vx) + y2 *(vy)) 
                            *exp(-width/(vx^2 + vy^2) * ( (x2 *(vy) - y2 *(vx))^2 + (x2* (vx) + y2 *(vy))^2 * γ^2 ))    )
        end
    end
end




#omegaIC calculates the Ω^2 matrix used in CQC calculations.
#And it uses to return sqrt(Ω) and inverse of sqrt(Ω) matrices to be used in the initial conditions
export omegaIC
function omegaIC(ϕ_s,ψ_s)
    
    #Calculate the matrix Ω^2
    Ω=zeros(N^2,N^2)
    # #! OMEGA TESTER -- Keep it for now in case need to test again.
        # CCD=zeros(N^2,N^2)
        # Ã=zeros(N^2,N^2)
        # B̃=zeros(N^2,N^2)
        # A=zeros(N^2,N^2)
        # B=zeros(N^2,N^2)
        # W=zeros(N^2,N^2)
        # Q=zeros(N^2,N^2)

    for J=1:N^2
        for K=1:N^2
            if (J==K)   #C, C̃, D
                Ω[J,K] = Ω[J,K]+ 2/dx^2 + 2/dy^2 + (m_ρ^2 + α*abs2(ϕ_s[J]) + β*ψ_s[K]^2)
                # CCD[J,K]= CCD[J,K]+ 2/dx^2 + 2/dy^2 + (m_ρ^2 + alpha*abs2(ϕ_s[J]) + beta*ψ_s[K]^2) #!

            elseif (J==K+1)         #Ã  and W
                Ω[J,K] = Ω[J,K] + -1/dy^2   #Ã
                # Ã[J,K] = Ã[J,K] + -1/dy^2 #!

                if (mod(K,N)==0 && K!=N^2)  #W
                    Ω[J,K] = Ω[J,K] + 1/dy^2
                    # W[J,K] = W[J,K] + 1/dy^2#!
                end

            elseif (J==1 && K==N^2) #Ã  and W
                Ω[J,K] = Ω[J,K] + -1/dy^2   #Ã 
                # Ã[J,K] = Ã[J,K] + -1/dy^2#!

                Ω[J,K] = Ω[J,K] + 1/dy^2     #W
                # W[J,K] = W[J,K] + 1/dy^2#!

            elseif (J==K-1)         #B̃  and W
                Ω[J,K] = Ω[J,K] + -1/dy^2   #B̃
                # B̃[J,K] = B̃[J,K] + -1/dy^2#!

                if (mod(J,N)==0 &&  J!=N^2) #W
                    Ω[J,K] = Ω[J,K] + 1/dy^2
                    # W[J,K] = W[J,K] + 1/dy^2#!
                end

            elseif (J==N^2 && K==1) #B̃ and W
                Ω[J,K] = Ω[J,K] + -1/dy^2   #B̃
                # B̃[J,K] = B̃[J,K] + -1/dy^2#!
                
                Ω[J,K] = Ω[J,K] + 1/dy^2 #W
                # W[J,K] = W[J,K] + 1/dy^2#!

            elseif (J==K+N && J>N)  #A
                Ω[J,K] = Ω[J,K] + -1/dx^2
                # A[J,K] = A[J,K] + -1/dx^2#!
            elseif (J in 1:N && K == J+N^2-N )  #A
                    Ω[J,K] = Ω[J,K] + -1/dx^2
                    # A[J,K] = A[J,K] + -1/dx^2#!

            elseif (J==K-N && J<N^2-N+1)     #B
                Ω[J,K] = Ω[J,K] + -1/dx^2
                # B[J,K] = B[J,K] + -1/dx^2 #!
            elseif (J in N^2-N+1:N^2 && K==J-(N^2-N))   #B
                    Ω[J,K] = Ω[J,K] + -1/dx^2
                    # B[J,K] = B[J,K] + -1/dx^2#!

            elseif (K==J-(N-1) && mod(J,N)==0)  #Q
                Ω[J,K] = Ω[J,K] + -1/dy^2
                # Q[J,K] = Q[J,K] + -1/dy^2#!
            elseif (J==K-(N-1) && mod(K,N)==0)  #Q
                Ω[J,K] = Ω[J,K] + -1/dy^2
                # Q[J,K] = Q[J,K] + -1/dy^2#!
            end
        end
    end

    ###################--CHOOSE ONE OF THEM--###########################
    #Taking sqrt of matrix Ω (version 1) #!--->This looks faster for N=3 test (but for N=60,70 got slower and uses more ram)
    # eigenValues=eigvals(Ω)
    # similarityMatrix = eigvecs(Ω)
    # diag_S_Ωzero=zeros(1:N^2,1:N^2)
    # for J=1:N^2
    #     for K=1:N^2
    #         if J==K
    #             diag_S_Ωzero[J,K] = eigenValues[J]^(1/4)
    #         end
    #     end
    # end
    # S_Ωzero= similarityMatrix * diag_S_Ωzero * transpose(similarityMatrix)
    # inv_S_Ωzero = inv(S_Ωzero)


    #Version 2 #!---> seems to be slower at the N=3 test (but for N=60,70 got faster and uses less ram!)
    S_Ωzero = sqrt(sqrt(Ω))  #S_ stands for square root
    inv_S_Ωzero = inv(S_Ωzero)

    ####################################################################



    return S_Ωzero, inv_S_Ωzero
    #!
    # return Ω,CCD,Ã,B̃,A,B,W,Q

end




export renormalization
    #Renormalization is done using the <ρ^2>_0 factor which we calculate here.
    #<ρ^2>_0 is calculated using Z values when |ϕ|=η and ψ=0. 
    #Here we use the shortcut that we assume at the boundaries of our lattice,
    #the fields are at their vacuum (|ϕ|=η and ψ=0) to begin with.
    #Thus, instead of calculating Z at each lattice point for vacuum values of ϕ and ψ,
    #we only use the value of Z from the boundary.
    #More specifically, bottom right corner of lattice is used (could be any point on boundary)
function renormalization(Z)

    #Calculate on the corner chunk
    if myrank==nprocs-1
        meanSqrRenorm = sum(abs2,Z[rx_l,:,ry_l,:])
    else
        meanSqrRenorm = 0.0
    end

    #Send to all
    meanSqrRenorm = MPI.Bcast(meanSqrRenorm, nprocs-1, comm) 

    if myrank==1
        @show meanSqrRenorm
    end
    
return meanSqrRenorm
end




export zeroPointEnergy
function zeroPointEnergy(Z,dZdt)

    if myrank==nprocs-1
        kEZRen = 0.
        gEZRen = 0.
        pEZRen = 0.
        for l=1:Nx
            for m=1:Ny
                kEZRen = kEZRen + ( abs2( dZdt[rx_l,l,ry_l,m,1] ) )/2
                gEZRen = gEZRen + ( abs2( (Z[rx_l,l,ry_l,m] - Z[rx_l-1,l,ry_l,m])/dx ) 
                                  + abs2( (Z[rx_l,l,ry_l,m] - Z[rx_l,l,ry_l-1,m])/dy ) )/2
                pEZRen = pEZRen + m_ρ^2*( abs2(Z[rx_l,l,ry_l,m,1]) )/2
            end
        end
        zPE = (kEZRen+gEZRen+pEZRen)/(dx*dy)
    else 
        zPE = 0.0
    end
    zPE = MPI.Bcast(zPE, nprocs-1, comm)

    if myrank==0
        @show zPE
    end

return zPE
end

end #module
module IC

include("Parameters.jl")
using .Parameters
using LinearAlgebra
include("IndexMap.jl")
using .IndexMap

#------------------------------------------------------------------------------------------------#
    #initialConditions sets the initial conditions for the fields ϕ,ψ,Z and calculates the renormalization factor
export initialConditions!
#Set initial conditions for the fields ϕ, ψ, ρ
function initialConditions!(ϕ,ψ,Z,dϕdt,dψdt,dZdt)
width=0.5
amp=10
    #ϕ and ψ i.c.
    for j=lx:rx
        x=j*dx
        for k=ly:ry
            y=k*dy
            ϕ[j,k,0] = eta
            ψ[j,k,0] = amp*exp(-x^2*width)*exp(-y^2*width)
            dϕdt[j,k,0] = 0
            dψdt[j,k,0] = 0
        end
    end

    #2-index to 1-index mapping
    ϕ_s = flattenDimension(ϕ[:,:,0])
    ψ_s = flattenDimension(ψ[:,:,0])
    #Get Sqrt(Omega) and its inverse matrices
    S_Ωzero, inv_S_Ωzero = omegaIC(ϕ_s,ψ_s)

    #Z (ρ) i.c.
    for J=1:N^2
        for K=1:N^2
            Z[J,K,1]=-im/sqrt(2) * inv_S_Ωzero[J,K]
            dZdt[J,K,1] = 1/sqrt(2) * S_Ωzero[J,K]
        end
    end

end
#------------------------------------------------------------------------------------------------#


#------------------------------------------------------------------------------------------------#
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
#------------------------------------------------------------------------------------------------#


#------------------------------------------------------------------------------------------------#
export renormalization
    #Renormalization is done using the <ρ^2>_0 factor which we calculate here.
    #<ρ^2>_0 is calculated using Z values when |ϕ|=η and ψ=0. 
    #Here we use the shortcut that we assume at the boundaries of our lattice,
    #the fields are at their vacuum (|ϕ|=η and ψ=0) to begin with.
    #Thus, instead of calculating Z at each lattice point for vacuum values of ϕ and ψ,
    #we only use the value of Z from the boundary.
    #More specifically, bottom right corner of lattice is used (could be any point on boundary)
function renormalization(Z)
    meanSqrRenorm = 0
    for K=1:N^2
        meanSqrRenorm = meanSqrRenorm + abs2(Z[N^2,K,1]) #There is an overall factor of 1/(dx*dy) which we omit here;
    end                                                  #It is added wherever we use this two-point function

    return meanSqrRenorm
end
#------------------------------------------------------------------------------------------------#

end
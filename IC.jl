module IC

include("Parameters.jl")
using .Parameters
using LinearAlgebra
include("IndexMap.jl")
using .IndexMap

#!
using OffsetArrays
#!
using Profile
using PProf

#------------------------------------------------------------------------------------------------#
    #initialConditions sets the initial conditions for the fields ϕ,ψ,Z and calculates the renormalization factor
export initialConditions!
#Set initial conditions for the fields ϕ, ψ, ρ
function initialConditions!(ϕ,ψ,Z,dϕdt,dψdt,dZdt)

    width=2.0
    amp=10
    # vel=0.5
    vx=0.3
    vy=0.4
    r0=1.25
    γ=1/sqrt(1-(vx^2+vy^2))
        #ϕ and ψ i.c.
        for j in axes(ϕ,1)
            x=j*dx
            x1 = x-r0
            x2 = x+r0
            for k in axes(ϕ,2)
                y=k*dy
                y1 = y-r0
                y2 = y+r0
                ϕ[j,k,0] = η 
                dϕdt[j,k,0] = 0.
                ψ[j,k,0] = amp*(exp( -width/(vx^2 + vy^2)
                                    * ( (x1 *(-vy) - y1 *(-vx))^2 + (x1* (-vx) + y1*(-vy))^2 * γ^2 ) )
                            + exp( -width/(vx^2 + vy^2) 
                                    * ( (x2 *vy - y2 *vx)^2 + (x2* vx + y2 *vy)^2 * γ^2 )) )
                dψdt[j,k,0] = ( 2amp *width *γ^2 *(x1 *(-vx) + y1 *(-vy)) 
                                *exp(-width/(vx^2 + vy^2) * ( (x1 *(-vy) - y1 *(-vx))^2 + (x1* (-vx) + y1 *(-vy))^2 * γ^2 ))
                              + 2amp *width *γ^2 *(x2 *(vx) + y2 *(vy)) 
                                *exp(-width/(vx^2 + vy^2) * ( (x2 *(vy) - y2 *(vx))^2 + (x2* (vx) + y2 *(vy))^2 * γ^2 ))    )
            end
        end

        #2-index to 1-index mapping
        ϕ_s = @views flattenDimension(ϕ[:,:,0])
        ψ_s = @views flattenDimension(ψ[:,:,0])
        #Get Sqrt(Omega) and its inverse matrices
        S_Ωzero, inv_S_Ωzero = omegaIC(ϕ_s,ψ_s)

        #Z (ρ) i.c.
        for J=1:N^2
            for K=1:N^2
                Z[J,K,1]=-im/sqrt(2) * inv_S_Ωzero[J,K]
                dZdt[J,K,1] = 1/sqrt(2) * S_Ωzero[J,K]
            end
        end

        #4-index notation (for omega as well)
        # #!
        # #Test 4-indexed S_Ωzero, inv_S_Ωzero
        # S_Ωzero_f = mapZTo4Index(S_Ωzero)
        # inv_S_Ωzero_f =mapZTo4Index(inv_S_Ωzero)
        #     #Z (ρ) i.c.
        #     @views Z[:,:,:,:,0] .= -im/sqrt(2) .* inv_S_Ωzero_f
        #     @views dZdt[:,:,:,:,0].= 1/sqrt(2) .* S_Ωzero_f
        # #!




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
    # # Taking sqrt of matrix Ω (version 1) #!--->This looks faster for N=3 test (but for N=60,70 got slower and uses more ram)
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

    meanSqrRenorm_v2 = sum(abs2,Z[N^2,:,1])#!
    println(meanSqrRenorm)#!
    println(meanSqrRenorm_v2)#!
    return meanSqrRenorm
end
#------------------------------------------------------------------------------------------------#

export zeroPointEnergy
function zeroPointEnergy(Z,dZdt)
    kEZRen = 0.
    gEZRen = 0.
    pEZRen = 0.
    for K=1:N^2
        kEZRen = kEZRen + ( abs2( dZdt[N^2,K,1] ) )/2
        gEZRen = gEZRen + ( abs2( (Z[N^2,K,1] - Z[N^2-N,K,1])/dx ) 
                          + abs2( (Z[N^2,K,1] - Z[N^2-1,K,1])/dy ) )/2
        pEZRen = pEZRen + m_ρ^2*( abs2(Z[N^2,K,1]) )/2
        # #!
        # kEZRen = kEZRen + ( abs2( dZdtREN[K] ) )/2
        # gEZRen = gEZRen + ( abs2( (ZREN[K,1] - ZREN[K,2])/dx ) 
        #                     + abs2( (ZREN[K,1] - ZREN[K,3])/dy ) )/2
        # pEZRen = pEZRen + m_ρ^2*( abs2(ZREN[K,1]) )/2
        # #!
    end

    zPE = (kEZRen+gEZRen+pEZRen)/(dx*dy)

    return zPE
end

end
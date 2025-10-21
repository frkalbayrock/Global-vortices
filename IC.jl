module IC

include("Parameters.jl")
using .Parameters
using LinearAlgebra
include("IndexMap.jl")
using .IndexMap
include("Auxiliary.jl")
using .Auxiliary_Routines
using Distributed
include("lapack_wrappers.jl")



#------------------------------------------------------------------------------------------------#
#initialConditions sets the initial conditions for the fields ϕ,ψ,Z and calculates the renormalization factor
export initialConditions!
#Set initial conditions for the fields ϕ, ψ, ρ
function initialConditions!(ϕ_gl,ψ_gl)

    #---Find the separation between Gaussians
    # r0 = gaussian_separation()
    r0 = 1.25#!

    #---Set i.c. for ϕ and ψ locally
    @everywhere workers() ic_ϕ_ψ!(ϕ,ψ,dϕdt,dψdt,$r0)


    #---Collect ϕ and ψ to global fields for Ω calculation
    for i=2:nprocs()
        lx_p, rx_p, ly_p, ry_p = chunker(i) #_p: physical
        ϕ_gl[lx_p:rx_p,ly_p:ry_p] .= (@fetchfrom i Main.ϕ[padd+1:Nx_loc+padd,
                                                            1+padd:padd+Ny_loc])::Array{ComplexF64, 2} 
        ψ_gl[lx_p:rx_p,ly_p:ry_p] .= (@fetchfrom i Main.ψ[padd+1:Nx_loc+padd,
                                                            1+padd:padd+Ny_loc])::Array{Float64, 2}        
    end


    #---Calculate Ω
    #2-index to 1-index mapping
    ϕ_s = flattenDimension(ϕ_gl)
    ψ_s = flattenDimension(ψ_gl)
    #Get Sqrt(Omega) and its inverse matrices
    S_Ωzero, inv_S_Ωzero = omegaIC(ϕ_s,ψ_s)


    #---Set i.c. for global Z
    #4-indexed Z and Ω
    S_Ωzero_f = mapZTo4Index(S_Ωzero)
    inv_S_Ωzero_f = mapZTo4Index(inv_S_Ωzero)



    #---Chunk Z and send it to workers
    @sync for p in workers()
        #Get the slices from sqrt(Ω) and sqrt(Ω)' to be sent to worker p
        @inline slices = get_worker_slices(p,inv_S_Ωzero_f,S_Ωzero_f)
        @spawnat p Main.send_and_update(slices)
    end

end






function gaussian_separation() 

    #Create the guide Gaussian for determining separation
    #from the given set of parameters for the wavepackets
    # ψGuide = OffsetArray(zeros(Float64, Nx,Ny),lx:rx,ly:ry)
    ψGuide = zeros(Float64, Nx,Ny)
    for j in axes(ψGuide,1)
        x = (j - Int(Nx/2-1)) * dx
        for k in axes(ψGuide,2)
            y = (k - Int(Nx/2-1)) * dy
            ψGuide[j,k] = amp*(exp( -width/(v^2)
                                * ( ( - y *(-v))^2 + (x* (-v) )^2 * γ^2 ) ))
        end
    end


    #Determine threshold from the amplitude
    threshold = 10^-3 * amp 

    #Find the position of the peak
    peak = argmax(ψGuide)

    #Find the position of cut-off
    cutoff = findfirst(y -> y > threshold, ψGuide[:, peak[2]])

    #Find the distance from cut-off to peak
    distance = dx * abs(cutoff - peak[1])

    #Determine r0
    r0 = 1/sqrt(2) * (distance + dx/2)
    @show r0

return r0
end




@everywhere workers() function send_and_update(data)
    Z[padd+1:Nx_loc+padd, :, 
        padd+1:Ny_loc+padd, :] .= -im/sqrt(2) .* data.inv_s_slice
    dZdt[:,:,:,:,1]            .=   1/sqrt(2) .* data.s_slice
end



@inline function get_worker_slices(p,inv_S_Ωzero_f,S_Ωzero_f)
    #Get the global ends of the lattice chunk
        #chunker() is written with offset arrays in mind.
        #Thus we shift the coordinates for base-1 arrays 
        #(since these sqrt and inverse sqrt matrices are not offset unlike fields)
    lx_p, rx_p, ly_p, ry_p = chunker(p)
    lx_p = Int(Nx/2+lx_p)
    rx_p = Int(Nx/2+rx_p)
    ly_p = Int(Ny/2+ly_p)
    ry_p = Int(Ny/2+ry_p)

    return (
        inv_s_slice = inv_S_Ωzero_f[lx_p:rx_p, :, ly_p:ry_p, :],
        s_slice     =     S_Ωzero_f[lx_p:rx_p, :, ly_p:ry_p, :]
    )
end




@everywhere workers() function ic_ϕ_ψ!(ϕ,ψ,dϕdt,dψdt,r0)

    #Find the chunk's physical coordinates and physical ends
    lx_p, rx_p, ly_p, ry_p = chunker(myid())


    #ϕ fluctuation parameters
    κ=2*2pi/L
    f_amp = 0.25

    #ϕ and ψ i.c.
    for j in padd+1:Nx_loc+padd
        x = (lx_p+(j-padd)-1)*dx
        x1 = x4 = x-r0
        x2 = x3 = x+r0
        for k=padd+1:Ny_loc+padd
            y = (ly_p+(k-padd)-1)*dy
            y1 = y3 = y-r0
            y2 = y4 = y+r0
            δϕ = f_amp*sin(κ*x)sin(κ*y) 
            ϕ[j,k] = η #+ im*δϕ
            dϕdt[j-padd,k-padd,1] = 0
            # ψ[j,k] = amp*(exp( -width/(vx^2 + vy^2)
            #                    * ( (x1 *(-vy) - y1 *(-vx))^2 + (x1* (-vx) + y1*(-vy))^2 * γ^2 ) )       
            #               + exp( -width/(vx^2 + vy^2) 
            #                  * ( (x2 *vy - y2 *vx)^2 + (x2* vx + y2 *vy)^2 * γ^2 )) )
            # dψdt[j-padd,k-padd,1] = ( 2amp *width *γ^2 *(x1 *(-vx) + y1 *(-vy)) 
            #                 *exp(-width/(vx^2 + vy^2) * ( (x1 *(-vy) - y1 *(-vx))^2 + (x1* (-vx) + y1 *(-vy))^2 * γ^2 ))
            #                 + 2amp *width *γ^2 *(x2 *(vx) + y2 *(vy)) 
            #                 *exp(-width/(vx^2 + vy^2) * ( (x2 *(vy) - y2 *(vx))^2 + (x2* (vx) + y2 *(vy))^2 * γ^2 ))    )
            ψ[j,k] = amp*(
                          exp( -width/(vx^2 + vy^2)
                                * ( (x1 *(-vy) - y1 *(-vx))^2 + (x1* (-vx) + y1*(-vy))^2 * γ^2 ) )  #-1-     
                        + exp( -width/(vx^2 + vy^2)     
                                * ( (x2 *(+vy) - y2 *(+vx))^2 + (x2* (+vx) + y2 *(+vy))^2 * γ^2 ))  #-2-
                        # + exp( -width/(vx^2 + vy^2)
                        #         * ( (x3 *(-vy) - y3 *(+vx))^2 + (x3* (+vx) + y3*(-vy))^2 * γ^2 ) )  #-3-
                        # + exp( -width/(vx^2 + vy^2)
                        #         * ( (x4 *(+vy) - y4 *(-vx))^2 + (x4* (-vx) + y4*(+vy))^2 * γ^2 ) )  #-4-
                        )

            dψdt[j-padd,k-padd,1] = 2amp *width *γ^2 *( 
                                      (x1 *(-vx) + y1 *(-vy)) 
                                        * exp(-width/(vx^2 + vy^2) * ( (x1 *(-vy) - y1 *(-vx))^2 + (x1* (-vx) + y1 *(-vy))^2 * γ^2 )) #-1-
                                    + (x2 *(+vx) + y2 *(+vy)) 
                                        * exp(-width/(vx^2 + vy^2) * ( (x2 *(+vy) - y2 *(+vx))^2 + (x2* (+vx) + y2 *(+vy))^2 * γ^2 )) #-2-
                                    # + (x3 *(+vx) + y3 *(-vy)) 
                                    #     * exp(-width/(vx^2 + vy^2) * ( (x3 *(-vy) - y3 *(+vx))^2 + (x3* (+vx) + y3 *(-vy))^2 * γ^2 )) #-3-
                                    # + (x4 *(-vx) + y4 *(+vy)) 
                                    #     * exp(-width/(vx^2 + vy^2) * ( (x4 *(+vy) - y4 *(-vx))^2 + (x4* (-vx) + y4 *(+vy))^2 * γ^2 )) #-4-
                                    )

        end
    end
end



   
#omegaIC calculates the Ω^2 matrix used in CQC calculations.
#And it uses to return sqrt(Ω) and inverse of sqrt(Ω) matrices to be used in the initial conditions
function omegaIC(ϕ_s,ψ_s)
    
    #Calculate the matrix Ω^2

    #Initialize Ω^2
    Ω = zeros(Float64, N2, N2)

    #Assign the values based on CQC
    for J=1:N2
        for K=1:N2
            if (J==K)   #C, C̃, D
                Ω[J,K] = Ω[J,K]+ 2/dx^2 + 2/dy^2 + (m_ρ^2 + α*abs2(ϕ_s[J]) + β*ψ_s[K]^2)

            elseif (J==K+1)         #Ã  and W
                Ω[J,K] = Ω[J,K] + -1/dy^2   #Ã

                if (mod(K,N)==0 && K!=N2)  #W
                    Ω[J,K] = Ω[J,K] + 1/dy^2
                end

            elseif (J==1 && K==N2) #Ã  and W
                Ω[J,K] = Ω[J,K] + -1/dy^2   #Ã 

                Ω[J,K] = Ω[J,K] + 1/dy^2     #W

            elseif (J==K-1)         #B̃  and W
                Ω[J,K] = Ω[J,K] + -1/dy^2   #B̃

                if (mod(J,N)==0 &&  J!=N2) #W
                    Ω[J,K] = Ω[J,K] + 1/dy^2
                end

            elseif (J==N2 && K==1) #B̃ and W
                Ω[J,K] = Ω[J,K] + -1/dy^2   #B̃
                
                Ω[J,K] = Ω[J,K] + 1/dy^2 #W

            elseif (J==K+N && J>N)  #A
                Ω[J,K] = Ω[J,K] + -1/dx^2

            elseif (J in 1:N && K == J+N2-N )  #A
                    Ω[J,K] = Ω[J,K] + -1/dx^2

            elseif (J==K-N && J<N2-N+1)     #B
                Ω[J,K] = Ω[J,K] + -1/dx^2

            elseif (J in N2-N+1:N2 && K==J-(N2-N))   #B
                    Ω[J,K] = Ω[J,K] + -1/dx^2

            elseif (K==J-(N-1) && mod(J,N)==0)  #Q
                Ω[J,K] = Ω[J,K] + -1/dy^2

            elseif (J==K-(N-1) && mod(K,N)==0)  #Q
                Ω[J,K] = Ω[J,K] + -1/dy^2

            end
        end
    end

    #Letting Julia know this is a Symmetric type for performance
    # Ω = Symmetric(Ω)

    # ###################--CHOOSE ONE OF THEM--###########################
    # Taking sqrt of matrix Ω (version 1) #!--->This looks faster for N=3 test (but for N=60,70 got slower and uses more ram)
    # println("---BEGIN---")
    # @time begin
    #     eigenVals, eigenVecs = eigen!(Ω; sortby=nothing)
    #     S_Ωzero = Symmetric(eigenVecs * Diagonal(eigenVals.^(0.25)) * eigenVecs')
    #     inv_S_Ωzero = Symmetric(eigenVecs * Diagonal(eigenVals.^(-0.25)) * eigenVecs')
    # end
    # println("---END---")
  

    #Version 2 #!---> seems to be slower at the N=3 test (but for N=60,70 got faster and uses less ram!)
    # S_Ωzero = zeros(N2,N2)
    # inv_S_Ωzero = zeros(N2,N2)
    # S_Ωzero = sqrt(sqrt(Ω))  #S_ stands for square root
    # inv_S_Ωzero = inv(S_Ωzero)

    # Version 3 (fastest and the most robust/reliable so far)
    H, Q = sym_tridiagonalize!(Ω)
    eigenVals, Z_T = stedc_manual(H; compute_vectors=true)
    eigenVecs = Q * Z_T
    S_Ωzero     = Symmetric(eigenVecs * Diagonal(eigenVals.^(0.25)) * eigenVecs')
    inv_S_Ωzero = Symmetric(eigenVecs * Diagonal(eigenVals.^(-0.25)) * eigenVecs')
    ####################################################################

return S_Ωzero, inv_S_Ωzero
end




export renormalization
    #Renormalization is done using the <ρ^2>_0 factor which we calculate here.
    #<ρ^2>_0 is calculated using Z values when |ϕ|=η and ψ=0. 
    #Here we use the shortcut that we assume at the boundaries of our lattice,
    #the fields are at their vacuum (|ϕ|=η and ψ=0) to begin with.
    #Thus, instead of calculating Z at each lattice point for vacuum values of ϕ and ψ,
    #we only use the value of Z from the boundary.
    #More specifically, bottom right corner of lattice is used (could be any point on boundary)
function renormalization()
    corner_id = nprocs()
    # Int(nprocs_perdim[1]+1)
    meanSqrRenorm = (@fetchfrom corner_id sum(abs2,Main.Z[Nx_loc+padd,:,Ny_loc+padd,:]))::Float64
    println("<ρ^2>: ",meanSqrRenorm)

return meanSqrRenorm
end


export zeroPointEnergy
function zeroPointEnergy()

    #Get the proc id of the corner of x_max,y_max
    corner_id = nprocs()
    #Fetch from the corner chunk
    Z_partial = (@fetchfrom corner_id Main.Z[rx_l-1:rx_l,:,ry_l-1:ry_l,:])::Array{ComplexF64, 4}
    dZdt_partial = (@fetchfrom corner_id Main.dZdt[rx_l-padd,:,ry_l-padd,:,1])::Array{ComplexF64, 2}

    #Calculate zPE
    kEZRen = 0.
    gEZRen = 0.
    pEZRen = 0.
    for l=1:Nx
        for m=1:Ny
            kEZRen +=  ( abs2( dZdt_partial[l,m] ) )/2
            gEZRen +=  ( abs2( (Z_partial[2,l,2,m] - Z_partial[1,l,2,m])/dx ) 
                                + abs2( (Z_partial[2,l,2,m] - Z_partial[2,l,1,m])/dy ) )/2
            pEZRen +=  m_ρ^2*( abs2(Z_partial[2,l,2,m,1]) )/2
        end
    end
    zPE = (kEZRen+gEZRen+pEZRen)/(dx*dy)
    println("zPE= ",zPE)
return zPE
end

end #module
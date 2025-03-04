#Production of Vortices with Quantum Mediation

include("Parameters.jl")
using .Parameters
include("IC.jl")
using .IC
include("Auxiliary.jl")
using .Auxiliary_Routines
# include("energy.jl")
# using .energy
include("IndexMap.jl")
using .IndexMap
include("Evolution.jl")
using .Time_Evolution
include("Constraints.jl")
using .Constraints_Conserveds
using OffsetArrays
using DelimitedFiles
using Plots; pythonplot()
using Printf

#!
using BenchmarkTools
using Profile
using PProf



function run_ev()

    #--Initilize the Field Arrays and
    #Standards:  ϕ and ψ are in 2D lattice // Z is flattened 1D N^2 lattice (for now)
    ϕ = im*zeros(Nx,Ny,2)
    ψ = zeros(Nx,Ny,2)
    Z = Array{ComplexF64,3}(undef, Nx^2,Ny^2,2) #This way seems to be faster and memory friendely for very large complex arrays. 
    dϕdt = im*zeros(Nx,Ny,2)
    dψdt = zeros(Nx,Ny,2)
    dZdt = Array{ComplexF64,3}(undef, Nx^2,Ny^2,2) #This way seems to be faster and memory friendely for very large complex arrays.

    #-Offset arrays for the symmetric lattice coordinates (for more natural physical indexing)
    ϕ = OffsetArray(ϕ,lx:rx,ly:ry,0:1)
    ψ = OffsetArray(ψ,lx:rx,ly:ry,0:1)
    dϕdt = OffsetArray(dϕdt,lx:rx,ly:ry,0:1)
    dψdt = OffsetArray(dψdt,lx:rx,ly:ry,0:1)


    #--Information about the run
    open("data/info.dat","w") do io
        #--Write Info For Graphs--!
        write(io,N,"x",N,"\n")  #number of lattice points
        writedlm(io,dx)         #lattice spacing
        writedlm(io,dt)         #time spacing
        writedlm(io,nt)         #number of time steps
        writedlm(io,nsnaps)     #number of snapshots
    end
    #Info
    println("Number of lattice points: ",Nx," x ",Ny)
    println("Size L of lattice in x-direction: ", Nx*dx )
    println("Lattice Spacing dx and dy :",dx)
    println("Time Spacing dt: ",dt)
    println("Number of time steps: ",nt)
    println("Threads: ",Threads.nthreads())#!
    #Data files
    ioϕ=open("data/initial_phi.dat","w")
    ioψ=open("data/initial_psi.dat","w")



    #----Initial Conditions----#
    initialConditions!(ϕ,ψ,Z,dϕdt,dψdt,dZdt)

    #Record initial conditions
    writedlm(ioϕ,ϕ[:,:,0])
    writedlm(ioψ,ψ[:,:,0])

    #!
    # #Check constraints and conserved quantities
    # @views constraints_checker(Z[:,:,1],dZdt[:,:,1])
    # @views conserved_checker(Z[:,:,1],dZdt[:,:,1])

    #----Renormalization----#
    meanSqrRenorm = renormalization(Z)
    zPE = zeroPointEnergy(Z,dZdt)

    #----Initial Energy----#
    totalE, ZED = energy(ϕ,ψ,Z,dϕdt,dψdt,dZdt,meanSqrRenorm,zPE)
    println("Total initial energy: ", totalE)

    #----Time Evolution----#
    @time time_evolve!(ϕ,ψ,Z,dϕdt,dψdt,dZdt,meanSqrRenorm,zPE)





    #Close data files
    close(ioϕ)
    close(ioψ)

end

@time run_ev()


# end 
#END OF CODE





#!
# open("data/energies/energy-density.dat","w") do io
#     Edensity_two = ravelDimension(Edensity)
#     println(typeof(Edensity_two))
#     writedlm(io,Edensity_two)
# end

#!
# Z_f = Array{ComplexF64,4}(undef, Nx,Nx,Ny,Ny)
# Z_f = mapZTo4Index(Z[:,:,1]) 
# @time Z[:,:,1] .= mapZTo2Index(Z_f) 


#!
#1d lattice
# ZED = ravelDimension(ZED_s)
#2d lattice 
# ZED = energy(ϕ,ψ,Z,dϕdt,dψdt,dZdt,meanSqrRenorm)

# #!
# open("data/energies/Z-ED.dat","w") do io
#     writedlm(io,ZED)
# end
# open("data/psi.dat","w") do io
#     writedlm(io,ψ[:,:,0])
# end

#!
# # meanSqrRenorm = 18.33
# width=0.05
# amp=10
# vel=0.5
# vx=vel
# vy=vel
# γ=1/sqrt(1-vel^2)
# x0=5
# for j=lx:rx
#     x=j*dx
#     for k=lx:rx
#         y=k*dy
#         # #Boost in y-direction
#         # ψ[j,k,0] =  amp*exp(-(γ*y)^2*width)*exp(-x^2*width)
#         # dψdt[j,k,0] = 2amp *width *vel *γ^2 *y *exp(-(γ*y)^2*width)*exp(-x^2*width)
#         # #Boost in x-direction
#         # ψ[j,k,0] =  amp*exp(-(γ*x)^2*width)*exp(-y^2*width)
#         # dψdt[j,k,0] = 2amp *width *vel *γ^2 *x *exp(-(γ*x)^2*width)*exp(-y^2*width)
#         #Boost in x and y-direction
#         ψ[j,k,0] = amp*exp(-width/(vx^2 + vy^2) * ( (x *vy - y *vx)^2 + (x* vx + y *vy)^2 * γ^2 ))
#         dψdt[j,k,0] = ( 2amp *width *γ^2 *(x *vx + y*vy) 
#                         *exp(-width/(vx^2 + vy^2) * ( (x *vy - y *vx)^2 + (x* vx + y *vy)^2 * γ^2 ))  )
#     end
# end
# x = collect(lx*dx:dx:rx*dx)
# y = collect(ly*dy:dy:ry*dy)
# # psi2 = OffsetArray(ψ[:,:,0],1:Nx,1:Ny)
# psi2 = OffsetArray(dψdt[:,:,0],1:Nx,1:Ny)
# surface(x,y,psi2[:,:]',zlims=(-1,1),xlabel="x",ylabel="y",zlabel="ψ")
# # plot(x,psi2[:,25],xlabel="x",ylabel="ψ")
# gui()
# savefig("PLotting/plots/x-y-boosted-psi-dot.png")
# # @time begin
#     # @btime time_evolve!(ϕ,ψ,Z,dϕdt,dψdt,dZdt,meanSqrRenorm)
#     # io = open("data/psi.dat","w")
# # @time begin
# meanSqrRenorm = 18.33
# end
# end


# #!
# open("data/Z-1timestep.dat","w") do io
#     writedlm(io,Z[:,:,1])
#     writedlm(io,"\n")
#     writedlm(io,dZdt[:,:,1])
# end





# baban=zeros(Nx,Ny)
# x = y = collect(lx*dx:dx:rx*dx)
# surface(x, y, ZED, c=:viridis)
# surface!(x,y,baban,c=:viridis)
# gui()
# readline()
# savefig("plots/ZED.png")

#!-----TESTING-----
#!
# open("data/phi.dat","w") do io
#     writedlm(io,ϕ[:,:,0])
# end
# open("data/psi-2d.dat","w") do io2
#     writedlm(io2,ψ[:,:,0])
# end

# #!
# ψ_s = zeros(N^2,2)
# @time begin
#     ψ_s[:,1] =  @views flattenDimension(ψ[:,:,0])
# end
# ψ[:,:,0] .= ravelDimension(view(ψ_s,:,1))


# x = y = collect(lx*dx:dx:rx*dx)
# # surface(x, y, abs2.(ϕ[:,:,0]), c=:viridis)#, zlim=(-10,10))
# surface!(x,y,ψ[:,:,0]',c=:viridis, xlabel="x", ylabel="y")
# gui()
# readline()
# savefig("psi2.png")


# ϕ_s = zeros(N^2)
# ψ_s = zeros(N^2)

# Ωzero, omegaZero = omegaIC(ϕ_s,ψ_s)

# # open("data/sqrtTest.txt","w") do io
# #     # writedlm(io,round.(Int,Ωzero))
# #     writedlm(io,Ωzero)
# #     write(io,"\n")
# #     # writedlm(io,round.(Int,omegaZero))
# #     writedlm(io,omegaZero)
# # end

# #!OMEGA TESTER -- keep for now in case need to test again.
# # io=open("data/ananTest.txt","w")

# # writedlm("data/anan.txt",Ω)
# # writedlm(io,"CCD")
# # writedlm(io,round.(Int,CCD))
# # writedlm(io,"At")
# # writedlm(io,round.(Int,Ã))
# # writedlm(io,"Bt")
# # writedlm(io,round.(Int,B̃))
# # writedlm(io,"A")
# # writedlm(io,round.(Int,A))
# # writedlm(io,"B")
# # writedlm(io,round.(Int,B))
# # writedlm(io,"W")
# # writedlm(io,round.(Int,W))
# # writedlm(io,"Q")
# # writedlm(io,round.(Int,Q))

# # close(io)

# # open("data/omegaCompare.txt","w") do io
# #     write(io,"Ω\n")
# #     writedlm(io,round.(Int,Ω))
# #     write(io,"\nomegaSummed\n")
# #     omegaSummed = CCD.+Ã.+B̃.+A.+B.+W.+Q
# #     writedlm(io,round.(Int,omegaSummed))
# # end


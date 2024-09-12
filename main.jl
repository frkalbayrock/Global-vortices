#Production of Vortices with Quantum Mediation

include("Parameters.jl")
using .Parameters
using Distributed
addprocs(prod(nprocs_perdim),topology=:all_to_all,lazy=true)

include("IC.jl")
include("Energy.jl")
include("Auxiliary.jl")
include("IndexMap.jl")
include("Evolution.jl")
include("Constraints.jl")
using .IC
using .Energy
using .Auxiliary_Routines
using .IndexMap
using .Time_Evolution
using .Constraints_Conserveds
using OffsetArrays
using DelimitedFiles
using Plots; pythonplot()
using Printf
# !
# using BenchmarkTools
# using Profile
# using PProf



@everywhere workers() begin
    include("Parameters.jl")
    include("IC.jl")
    include("Energy.jl")
    include("Auxiliary.jl")
    include("IndexMap.jl")
    include("Evolution.jl")
    include("Constraints.jl")
    using .Parameters
    using .IC
    using .Energy
    using .Auxiliary_Routines
    using .IndexMap
    using .Time_Evolution
    using .Constraints_Conserveds
end

function run_ev()
    #--Initilize the Chunk Field Arrays and
    #Standards:  ϕ and ψ are in 2D lattice // Z can be on the flattened 1D N^2 lattice or native 2D lattice
    ϕ_gl =    im*zeros(Nx, Ny, 2)
    ψ_gl =       zeros(Nx, Ny, 2)
    ϕ_gl = OffsetArray(ϕ_gl,lx:rx,ly:ry,0:1)
    ψ_gl = OffsetArray(ψ_gl,lx:rx,ly:ry,0:1)
    # #2-index Z
    # # Z_gl =    Array{ComplexF64,3}(undef, Nx^2,Ny^2,2) #!WE PROBABLY DONT NEED TO DEFINE THESE HERE 
    # # dZdt_gl = Array{ComplexF64,3}(undef, Nx^2,Ny^2,2) #!NO, Z would be a problem. MAYBE NOT.......
    #4-index Z
    Z_gl =    Array{ComplexF64,5}(undef, Nx, Nx, Ny, Ny, 2)
    dZdt_gl = Array{ComplexF64,5}(undef, Nx, Nx, Ny, Ny, 2)

    #Define the field on chunks
    @sync @everywhere workers() begin 
        ϕ =    im*zeros(Nx_loc+padding_size, Ny_loc+padding_size, 2)
        ψ =       zeros(Nx_loc+padding_size, Ny_loc+padding_size, 2)
        dϕdt = im*zeros(Nx_loc+padding_size, Ny_loc+padding_size, 2) #!it looks like these don't need the padding.
        dψdt =    zeros(Nx_loc+padding_size, Ny_loc+padding_size, 2) #!but easier for "for" loops when ϕ and dϕdt are in the same one
        #2-index Z
        # Z =    Array{ComplexF64,3}(undef, Nx^2,Ny^2,2) #This way seems to be faster and memory friendely for very large complex arrays. 
        # dZdt = Array{ComplexF64,3}(undef, Nx^2,Ny^2,2)
        #4-index Z
        Z =    Array{ComplexF64,5}(undef, Nx_loc+padding_size, Nx, Ny_loc+padding_size, Ny, 2)
        dZdt = Array{ComplexF64,5}(undef, Nx_loc+padding_size, Nx, Ny_loc+padding_size, Ny, 2)
    end


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
    println("Procs: ",nprocs())
    #Data files
    ioϕ=open("data/initial_phi.dat","w")
    ioψ=open("data/initial_psi.dat","w")



    #----Initial Conditions----#
    @time initialConditions!(ϕ_gl,ψ_gl,Z_gl,dZdt_gl)

    #Record initial conditions
    writedlm(ioϕ,ϕ_gl[:,:,0])
    writedlm(ioψ,ψ_gl[:,:,0])


    # #!Turn on if you wanna check constraints - skipping for now
    # @views Z_t = mapZTo2Index(Z_gl[:,:,:,:,0])
    # @views dZdt_t = mapZTo2Index(dZdt_gl[:,:,:,:,0])
    # @views constraints_checker(Z_t[:,:,1],dZdt_t[:,:,1])
    # @views conserved_checker(Z_t[:,:,1],dZdt_t[:,:,1])


    #----Renormalization----#
    meanSqrRenorm = renormalization()
    zPE = zeroPointEnergy()


    #----Initial Energy----#
    totalE, ZED = energy(meanSqrRenorm,zPE)
    ZedIO = open("data/energies/ZED.dat","w")
    writedlm(ZedIO,ZED)
    println("Total initial energy: ",totalE)


    #----Time Evolution----#
    @time time_evolve!(ϕ_gl,ψ_gl,meanSqrRenorm,zPE)



    #Close data files
    close(ioϕ)
    close(ioψ)

end 
@time run_ev()
rmprocs(workers())
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


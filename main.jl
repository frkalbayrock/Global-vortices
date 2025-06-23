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
include("FindVortex.jl")
using .IC
using .Energy
using .Auxiliary_Routines
using .IndexMap
using .Time_Evolution
using .Constraints_Conserveds
using .FindVortex
using OffsetArrays
using DelimitedFiles
using Plots; pythonplot()
using Printf
# !
using BenchmarkTools
using Profile
using PProf
using LinearAlgebra



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


#Open the directories
run(`mkdir -p data`)
run(`mkdir -p data/energies`)



#---Initialize the field arrays

    #Initilize the Global Lattice Arrays
    b = @allocated begin #!we keep these now until we finish the testing
    const ϕ_gl = OffsetArray(im*zeros(Nx, Ny),lx:rx,ly:ry)
    const ψ_gl = OffsetArray(zeros(Nx, Ny),lx:rx,ly:ry)
    const ZED_gl = OffsetArray(zeros(Nx,Ny),lx:rx,ly:ry)
    end
    println("Allocated: ",b/1e6," MB")

    #Initilize the Local Chunk Field Arrays
        #Standards:  ϕ and ψ are in 2D lattice // Z can be on the flattened 1D N^2 lattice or native 2D lattice
    @everywhere workers() begin
        const ϕ =    im*zeros(Nx_loc+padding_size, Ny_loc+padding_size)
        const ψ =       zeros(Nx_loc+padding_size, Ny_loc+padding_size)
        const dϕdt = im*zeros(Nx_loc, Ny_loc, 2)
        const dψdt =    zeros(Nx_loc, Ny_loc, 2)
        const ZED  =    zeros(Nx_loc, Ny_loc)
        #2-index Z
        # Z =    Array{ComplexF64,3}(undef, Nx^2,Ny^2,2) #This way seems to be faster and memory friendely for very large complex arrays. 
        # dZdt = Array{ComplexF64,3}(undef, Nx^2,Ny^2,2)
        #4-index Z
        const Z =    Array{ComplexF64,4}(undef, Nx_loc+padding_size, Nx, Ny_loc+padding_size, Ny)
        const dZdt = Array{ComplexF64,5}(undef, Nx_loc, Nx, Ny_loc, Ny, 2)
    end



#---Find neighbors of the chunks in the domain decomposition topology
    @everywhere workers() begin
        const n_left, n_right, n_bottom, n_top = find_neighbours(myid())
    end


#---BLAS Multi-threading Control
    BLAS.set_num_threads(prod(nprocs_perdim))






#---Main run function
function run_ev()

    #--Information about the run
    open("data/info.dat","w") do io
        #--Write Info For Graphs--!
        println(io,Nx,"x",Ny)  #number of lattice points
        println(io,dx)         #lattice spacing
        println(io,dt)         #time spacing
        println(io,nt)         #number of time steps
        println(io,nsnaps)     #number of snapshots
    end
    #Info
    println("Number of lattice points: ",Nx," x ",Ny)
    println("Size L of lattice in x-direction: ", L )
    println("Lattice Spacing dx and dy :",dx)
    println("Time Spacing dt: ",dt)
    println("Number of time steps: ",nt)
    println("Threads: ",Threads.nthreads())#!
    println("BLAS-Threads: ",BLAS.get_num_threads())
    println("Procs: ",nprocs())
    #Initial Data files
    ioϕ=open("data/initial_phi.dat","w")
    ioψ=open("data/initial_psi.dat","w")
    ioZed=open("data/initial_ZED.dat","w")



    #----Initial Conditions----#
    a = @allocated begin #!we keep these now until we finish the testing
    initialConditions!(ϕ_gl,ψ_gl)
    end
    println("Initial conditions allocated: ",a/1e6," MB")

    #Record initial conditions
    writedlm(ioϕ,ϕ_gl[:,:])
    writedlm(ioψ,ψ_gl[:,:])


    # #!Turn on if you wanna check constraints - skipping for now
    # @views Z_t = mapZTo2Index(Z_gl[:,:,:,:])
    # @views dZdt_t = mapZTo2Index(dZdt_gl[:,:,:,:])
    # @views constraints_checker(Z_t[:,:],dZdt_t[:,:])
    # @views conserved_checker(Z_t[:,:],dZdt_t[:,:])


    #----Renormalization----#
    meanSqrRenorm = renormalization()
    zPE = zeroPointEnergy()


    #----Initial Energy----#
    totalE = energy!(ZED_gl,meanSqrRenorm,zPE)
    writedlm(ioZed,ZED_gl)
    println("Total initial energy: ",totalE)


    #----Time Evolution----#
    @time time_evolve!(ϕ_gl,ψ_gl,ZED_gl,meanSqrRenorm,zPE)
    # @btime time_evolve!($ϕ_gl,$ψ_gl,$ZED_gl,$meanSqrRenorm,$zPE)
    # display(@benchmark time_evolve!($ϕ_gl,$ψ_gl,$ZED_gl,$meanSqrRenorm,$zPE) samples=10)


    #----Check for Vortices----#
    vortex_pos, anti_vortex_pos = vortex_finder(ϕ_gl)   #!this may not be necessary.
 
    #Close data files
    close(ioϕ)
    close(ioψ)
    close(ioZed)

end 

@time run_ev()
#!
# Profile.Allocs.@profile sample_rate=0.01 begin
#     run_ev()
# end
# PProf.Allocs.pprof(from_c=false)
#!



rmprocs(workers())
println("Removing workers done.")
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


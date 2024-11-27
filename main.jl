#Production of Vortices with Quantum Mediation
using MPI
MPI.Init()

include("Parameters.jl")
using .Parameters
include("IC.jl")
using .IC
include("Auxiliary.jl")
using .Auxiliary_Routines
include("Energy.jl")
using .Energy
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

include("MPIAux.jl")
using .MPIAux

#!
using BenchmarkTools
using Profile
using PProf

#------MPI PART------#
const comm = MPI.COMM_WORLD
const myrank = MPI.Comm_rank(comm)
const nprocs = MPI.Comm_size(comm)
const comm_cart = MPI.Cart_create(comm, nprocs_perdim; periodic=periods, reorder=false)
const coords_cart = MPI.Cart_coords(comm_cart, myrank)
#!We'll put all of these stuff into a initializer module where all MPI stuff will reside to be called whenever needed
#!IDK yet if it's gonna work.
if myrank==0
    const ϕ_gl = OffsetArray(im*zeros(Nx,Ny),lx:rx,ly:ry)
    const ψ_gl = OffsetArray(zeros(Nx,Ny),lx:rx,ly:ry)
    const ZED_gl  = OffsetArray(zeros(Nx,Ny),lx:rx,ly:ry)
end

function run_ev()

    #--Information about the run
    if myrank==0
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
        println("Size L of lattice in x-direction: ", Nx*dx )
        println("Lattice Spacing dx and dy :",dx)
        println("Time Spacing dt: ",dt)
        println("Number of time steps: ",nt)
        println("Threads: ",Threads.nthreads())#!
        #Data files
        ioϕ=open("data/initial_phi.dat","w")
        ioψ=open("data/initial_psi.dat","w")
    end

    #--Initilize the Field Arrays and
    #Standards:  ϕ and ψ are in 2D lattice // Z is flattened 1D N^2 lattice (for now)
    ϕ =   im.* zeros(Nx_loc+padding,Ny_loc+padding)
    ψ =        zeros(Nx_loc+padding,Ny_loc+padding)
    dϕdt = im.*zeros(Nx_loc+padding,Ny_loc+padding,2)
    dψdt =     zeros(Nx_loc+padding,Ny_loc+padding,2)
    ZED =      zeros(Nx_loc+padding,Ny_loc+padding,2)
    # 2-index Z
    # Z = Array{ComplexF64,3}(undef, (Nx_loc)*(Ny_loc),Nx*Ny,2)
    # dZdt = Array{ComplexF64,3}(undef, (Nx_loc)*(Ny_loc),Nx*Ny,2)
    # 4-index Z
    Z =    Array{ComplexF64,4}(undef, Nx_loc+padding, Nx, Ny_loc+padding, Ny)
    dZdt = Array{ComplexF64,5}(undef, Nx_loc+padding, Nx, Ny_loc+padding, Ny, 2)


    #----Initial Conditions----#
    @time initialConditions!(ϕ,ψ,Z,dϕdt,dψdt,dZdt)

    if myrank==0
        #Record initial conditions
        writedlm(ioϕ,ϕ_gl[:,:])
        writedlm(ioψ,ψ_gl[:,:])
    end


    #!
    # #Check constraints and conserved quantities
    # @views constraints_checker(Z[:,:,1],dZdt[:,:,1])
    # @views conserved_checker(Z[:,:,1],dZdt[:,:,1])

    #----Renormalization----#
    meanSqrRenorm = renormalization(Z)
    zPE = zeroPointEnergy(Z,dZdt)


    #----Initial Energy----#
    totalE = energy(ϕ,ψ,Z,dϕdt,dψdt,dZdt,ZED,meanSqrRenorm,zPE)
    if myrank==0
        ZedIO = open("data/energies/ZED.dat","w")
        writedlm(ZedIO,ZED_gl)
        println("Total initial energy: ",totalE)
    end
    

    #----Time Evolution----#
    @time time_evolve!(ϕ,ψ,Z,dϕdt,dψdt,dZdt,ZED,meanSqrRenorm,zPE)



    #Close data files
    if myrank==0
        close(ioϕ)
        close(ioψ)
    end
    MPI.Finalize()
end 

@time run_ev()

#END OF CODE
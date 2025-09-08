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
    const ϕ_gl = OffsetArray(zeros(ComplexF64, Nx, Ny),lx:rx,ly:ry)
    const ψ_gl = OffsetArray(zeros(Float64, Nx, Ny),lx:rx,ly:ry)
    const ZED_gl = OffsetArray(zeros(Float64, Nx,Ny),lx:rx,ly:ry)

    #Initilize the Local Chunk Field Arrays
        #Standards:  ϕ and ψ are in 2D lattice // Z can be on the flattened 1D N^2 lattice or native 2D lattice
    @everywhere workers() begin
        const ϕ =    zeros(ComplexF64, Nx_loc+padding_size, Ny_loc+padding_size)
        const ψ =    zeros(Nx_loc+padding_size, Ny_loc+padding_size)
        const dϕdt = zeros(ComplexF64, Nx_loc, Ny_loc, 2)
        const dψdt = zeros(Nx_loc, Ny_loc, 2)
        const ZED  = zeros(Nx_loc, Ny_loc)
        #2-index Z
        # Z =    Array{ComplexF64,3}(undef, Nx^2,Ny^2,2) #This way seems to be faster and memory friendely for very large complex arrays. 
        # dZdt = Array{ComplexF64,3}(undef, Nx^2,Ny^2,2)
        #4-index Z
        const Z =    Array{ComplexF64,4}(undef, Nx_loc+padding_size, Nx, Ny_loc+padding_size, Ny)
        const dZdt = Array{ComplexF64,5}(undef, Nx_loc, Nx, Ny_loc, Ny, 2)
        #Flux fields for vectorized time evolution
        const ϕ_flux = Array{Float64,2}(undef, Nx_loc, Ny_loc)
        const ψ_flux = Array{Float64,2}(undef, Nx_loc, Ny_loc)
        const Z_flux = Array{ComplexF64,4}(undef, Nx_loc, Nx, Ny_loc, Ny)
        const meanSqr_Rho = Array{Float64,2}(undef, Nx_loc, Ny_loc)
    end



#---Find neighbors of the chunks in the domain decomposition topology
    @everywhere workers() begin
        const n_left, n_right, n_bottom, n_top = find_neighbours(myid())
    end


#---Initialize the "Remote Channels" - networking between workers()
    initialize_channels()


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
    initialConditions!(ϕ_gl,ψ_gl)

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
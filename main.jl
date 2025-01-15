#Production of Vortices with Quantum Mediation

include("Parameters.jl")
using .Parameters
using Distributed
using DistributedArrays
addprocs(prod(nprocs_perdim),topology=:all_to_all,lazy=true)
@everywhere using DistributedArrays
@everywhere using DistributedArrays.SPMD

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
# using BenchmarkTools
using Profile
using PProf



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
    using DistributedArrays
end


#--Initilize the Global Lattice Arrays 
    #Standards:  ϕ and ψ are in 2D lattice // Z can be on the flattened 1D N^2 lattice or native 2D lattice
    const ϕ = DArray((Nx_padd,Ny_padd),workers(),[nprocs_perdim[1],nprocs_perdim[2]]) do I
        im*zeros(map(length,I)...)
    end
    const ψ = DArray((Nx_padd,Ny_padd),workers(),[nprocs_perdim[1],nprocs_perdim[2]]) do I
        zeros(map(length,I)...)
    end
    const dϕdt = DArray((Nx,Ny,2),workers(),[nprocs_perdim[1],nprocs_perdim[2],1]) do I
        im*zeros(map(length,I)...)
    end
    const dψdt = DArray((Nx,Ny,2),workers(),[nprocs_perdim[1],nprocs_perdim[2],1]) do I
           zeros(map(length,I)...)
    end
    const Z = DArray((Nx_padd,Nx,Ny_padd,Ny),workers(),[nprocs_perdim[1],1,nprocs_perdim[2],1]) do I
        im*zeros(length(I[1]),length(I[2]),length(I[3]),length(I[4]))
    end
    const dZdt = DArray((Nx,Nx,Ny,Ny,2),workers(),[nprocs_perdim[1],1,nprocs_perdim[2],1,1]) do I
        im*zeros(length(I[1]),length(I[2]),length(I[3]),length(I[4]),length(I[5]))
    end
    const ZED = DArray((Nx,Ny),workers(),[nprocs_perdim[1],nprocs_perdim[2]]) do I
        zeros(map(length,I)...)
    end
    #Global fields only for recording purposes --#!might change later if we can find an efficient way of doing it. 
    const ϕ_gl = im*zeros(Nx,Ny)
    const ψ_gl =    zeros(Nx,Ny)






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
    println("Procs: ",nprocs())
    #Data files
    ioϕ=open("data/initial_phi.dat","w")
    ioψ=open("data/initial_psi.dat","w")
    ioZED=open("data/initial_ZED.dat","w")


    #----Initial Conditions----#
    @time initialConditions!(ϕ,ψ,Z,dϕdt,dψdt,dZdt,ϕ_gl,ψ_gl)


    #Record initial conditions #!TEST IT FOR SPEED/ALLOCATIONS
    open("data/initial_phi.dat","w") do ioϕ
    open("data/initial_psi.dat","w") do ioψ
    for p in workers()
        lx_p, rx_p, ly_p, ry_p = distChunker(p)
        ϕ_gl[lx_p:rx_p,ly_p:ry_p] .= @fetchfrom p localpart(ϕ)[1+padd:Nx_loc+padd,1+padd:Ny_loc+padd,1]
        ψ_gl[lx_p:rx_p,ly_p:ry_p] .= @fetchfrom p localpart(ψ)[1+padd:Nx_loc+padd,1+padd:Ny_loc+padd,1]
    end
    writedlm(ioϕ,ϕ_gl[:,:])
    writedlm(ioψ,ψ_gl[:,:])
    # writedlm(ioψ,ψ[:,:,1]) #! this doesn't work for me since there are paddings 
    #                         #!either remove paddings on the plotting side or just use a separate "global" field array
    end
    end



    # #!Turn on if you wanna check constraints - skipping for now
    # @views Z_t = mapZTo2Index(Z_gl[:,:,:,:])
    # @views dZdt_t = mapZTo2Index(dZdt_gl[:,:,:,:])
    # @views constraints_checker(Z_t[:,:],dZdt_t[:,:])
    # @views conserved_checker(Z_t[:,:],dZdt_t[:,:])


    #----Renormalization----#
    meanSqrRenorm =  renormalization(Z)
    zPE = zeroPointEnergy(Z,dZdt)


    #----Initial Energy----#
    totalE = energy(ϕ,ψ,Z,dϕdt,dψdt,dZdt,ZED,meanSqrRenorm,zPE)
    writedlm(ioZED,ZED)
    println("Total initial energy: ",totalE)


    #----Time Evolution----#
    @time time_evolve!(ϕ,ψ,Z,dϕdt,dψdt,dZdt,ZED,meanSqrRenorm,zPE,ϕ_gl,ψ_gl)
    # Profile.Allocs.@profile sample_rate=0.01 time_evolve!(ϕ,ψ,Z,dϕdt,dψdt,dZdt,ZED,meanSqrRenorm,zPE,ϕ_gl,ψ_gl)
    # PProf.Allocs.pprof(from_c=false) 


end 

@time run_ev()
rmprocs(workers())
println("Removing workers done.")
#END OF CODE



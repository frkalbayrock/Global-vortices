#Production of Vortices with Quantum Mediation

include("Parameters.jl")
using .Parameters
using Distributed
addprocs(prod(nprocs_perdim),topology=:all_to_all,lazy=true)
using DistributedArrays

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
    const ϕ = im.*dzeros(Nx_padd,Ny_padd,2)
    const ψ =     dzeros(Nx_padd,Ny_padd,2)
    const dϕdt = im.*dzeros(Nx,Ny,2)
    const dψdt =     dzeros(Nx,Ny,2)
    const Z = DArray((Nx_padd,Nx,Ny_padd,Ny,2),workers(),[nprocs_perdim[1],1,nprocs_perdim[2],1,1]) do I
        im*zeros(length(I[1]),length(I[2]),length(I[3]),length(I[4]),length(I[5]))
    end
    const dZdt = DArray((Nx,Nx,Ny,Ny,2),workers(),[nprocs_perdim[1],1,nprocs_perdim[2],1,1]) do I
        im*zeros(length(I[1]),length(I[2]),length(I[3]),length(I[4]),length(I[5]))
    end
    const ZED = dzeros(Nx,Ny)
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
#    Profile.Allocs.@profile sample_rate=0.01 time_evolve!(ϕ,ψ,Z,dϕdt,dψdt,dZdt,ZED,meanSqrRenorm,zPE,ϕ_gl,ψ_gl)
#     PProf.Allocs.pprof(from_c=false) 


    #----Check for Vortices----#
    # vortex_pos, anti_vortex_pos = vortex_finder(ϕ)
 

    #Close data files
    close(ioϕ)
    close(ioψ)

end 

@time run_ev()
rmprocs(workers())
println("Removing workers done.")
#END OF CODE




# if myrank == 1
#     Profile.Allocs.@profile sample_rate=0.01 run_ev()
#     PProf.Allocs.pprof(from_c=false)
# else
#     @time run_ev()
# end



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


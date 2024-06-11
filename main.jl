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
using OffsetArrays
using DelimitedFiles
using Plots; pythonplot()
using Printf




#--Initilize the Field Arrays and
#Standards:  ϕ and ψ are in 2D lattice // Z is flattened 1D N^2 lattice (for now)
ϕ = im*zeros(Nx,Ny,2)
ψ = zeros(Nx,Ny,2)
Z = im*zeros(Nx^2,Nx^2,2)
dϕdt = im*zeros(Nx,Ny,2)
dψdt = zeros(Nx,Ny,2)
dZdt = im*zeros(Nx^2,Nx^2,2)

#-Offset arrays for the symmetric lattice coordinates (for more natural physical indexing)
ϕ = OffsetArray(ϕ,lx:rx,ly:ry,0:1)
ψ = OffsetArray(ψ,lx:rx,ly:ry,0:1)
dϕdt = OffsetArray(dϕdt,lx:rx,ly:ry,0:1)
dψdt = OffsetArray(dψdt,lx:rx,ly:ry,0:1)

#--Initial Conditions
initialConditions!(ϕ,ψ,Z,dϕdt,dψdt,dZdt)

#--Renormalization
meanSqrRenorm = renormalization(Z)

#--Initial Energy
ZED_s = energy(ϕ,ψ,Z,dϕdt,dψdt,dZdt,meanSqrRenorm)

#!
#1d lattice
ZED = ravelDimension(ZED_s)
# #2d lattice 
# ZED = energy(ϕ,ψ,Z,dϕdt,dψdt,dZdt,meanSqrRenorm)

#!
open("data/energies/Z-ED.dat","w") do io
    writedlm(io,ZED)
end
open("data/psi.dat","w") do io
    writedlm(io,ψ[:,:,0])
end


#--Time Evolution












# baban=zeros(Nx,Ny)
# x = y = collect(lx*dx:dx:rx*dx)
# surface(x, y, ZED, c=:viridis)
# surface!(x,y,baban,c=:viridis)
# gui()
# readline()
# savefig("plots/ZED.png")

# #!-----TESTING-----
# #!
# open("data/phi.dat","w") do io
#     writedlm(io,ϕ[:,:,0])
# end
# open("data/psi.dat","w") do io
#     writedlm(io,ψ[:,:,0])
# end
# x = y = collect(lx*dx:dx:rx*dx)
# # surface(x, y, abs2.(ϕ[:,:,0]), c=:viridis, zlim=(-10,10))
# surface!(x,y,ψ[:,:,0],c=:viridis)
# # gui()
# # readline()
# savefig("psi.png")


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


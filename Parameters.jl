module Parameters
export dx,dy,dt,N,Nx,Ny,L,lx,rx,ly,ry,nt,nsnaps, dt_half
export λ,α,β,η,m_ρ,m_ϕ,m_ψ
export ndims, nprocs_perdim, periods
export Nx_loc, Ny_loc, padding, padd
export lx_l,rx_l,ly_l,ry_l
export snapInterval

#--Lattice Parameters
const N=20
const Nx=N
const Ny=N
const lx=Int16.(-Nx/2+1)
const rx=Int16.(Nx/2)
const ly=lx
const ry=rx
const dx=0.4
const dy=dx
const dt=dx/50.0
const dt_half=dt/2
const L=Nx*dx
const nt=2000
#!maybe add L_x and L_y later (if needed)

#--Couplings and Masses #!TO BE FINE_TUNED
const λ=1
const α=0.5#1
const β=0.5#1
const η=1
const m_ρ=1
const m_ϕ=sqrt(λ*η^2)
const m_ψ=1

#--Data Recorder
const nsnaps=200
if round(nt/nsnaps,RoundDown) == 0
    const snapInterval = 1
else
    const snapInterval = round(Int,nt/nsnaps,RoundDown)
end

#--Parallel: MPI Parameters
const ndims = 2
const nprocs_perdim = [4 2]
const periods = [true true]
#--Parallel: Local array sizes and indices for OffsetArrays
const Nx_loc = Int(Nx/nprocs_perdim[1])
const Ny_loc = Int(Ny/nprocs_perdim[2])
#Padding to be added on top of the physical sizes
const padding = 2 #depends on the order of derivatives
const padd = Int(padding/2) #padding_perside
#Physically useful ends of local padded arrays
const lx_l = padd + 1
const rx_l = padd + Nx_loc
const ly_l = padd + 1
const ry_l = padd + Ny_loc



end #module
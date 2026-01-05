module Parameters
export dx,dy,dt,N,Nx,Ny,L,lx,rx,ly,ry,nt,nsnaps, dt_half,dx2,dy2,N2
export λ,α,β,η,m_ρ,m_ϕ,m_ψ
export nprocs_perdim,Nx_loc,Ny_loc,padding_size,padd
export lx_l,rx_l,ly_l,ry_l 
export width, amp, vx, vy, v, γ, vel

#--Lattice Parameters
const N=20
const Nx=N
const Ny=N
const N2=N^2
const lx=Int(-Nx/2+1)
const rx=Int(Nx/2)
const ly=lx
const ry=rx
const dx=0.4
const dy=dx
const dt=dx/50.0
const dt_half=dt/2
const dx2 = dx^2
const dy2 = dy^2
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

#--Gaussian Parameters
const width=2.0
const amp=20.0
const vel=0.4
const vx=vel/sqrt(2)
const vy=vel/sqrt(2)
const γ=1/sqrt(1-(vx^2+vy^2))

#--Data Recorder
const nsnaps=200

#--Parallel
const ndims = 2
const nprocs_perdim = [2 2]
#Physical sizes of chunks
const Nx_loc = Int(Nx/nprocs_perdim[1])
const Ny_loc = Int(Ny/nprocs_perdim[2])
#Padding to be added on top of the physical sizes
const padding_size = 2 #depends on the order of derivatives
const padd = Int(padding_size/2) #padding_perside
#Physically useful ends of local padded arrays
const lx_l = padd + 1
const rx_l = padd + Nx_loc
const ly_l = padd + 1
const ry_l = padd + Ny_loc


end #module
module Parameters
export dx,dy,dt,N,Nx,Ny,lx,rx,ly,ry,nt
export λ,α,β,eta,m_ρ,m_ϕ,m_ψ

#--Lattice Parameters
N=50
Nx=N
Ny=N
lx=Int16.(-Nx/2+1)
rx=Int16.(Nx/2)
ly=lx
ry=rx
dx=0.4
dy=dx
dt=dx/50.0
nt=1000
#!maybe add L_x and L_y later (if needed)

#--Couplings and Masses #!TO BE FINE_TUNED
λ=1
α=1
β=1
eta=1
m_ρ=1
m_ϕ=1
m_ψ=1

#--Data Recorder
nsnaps=20

end
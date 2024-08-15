module Parameters
export dx,dy,dt,N,Nx,Ny,lx,rx,ly,ry,nt,nsnaps
export λ,α,β,η,m_ρ,m_ϕ,m_ψ

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
const nt=2000
#!maybe add L_x and L_y later (if needed)

#--Couplings and Masses #!TO BE FINE_TUNED
const λ=1
const α=1
const β=1
const η=1
const m_ρ=1
const m_ϕ=sqrt(λ*η^2)
const m_ψ=1

#--Data Recorder
const nsnaps=200

# #--Lattice Parameters
# N=20
# Nx=N
# Ny=N
# lx=Int16.(-Nx/2+1)
# rx=Int16.(Nx/2)
# ly=lx
# ry=rx
# dx=0.4
# dy=dx
# dt=dx/50.0
# nt=1000
# #!maybe add L_x and L_y later (if needed)

# #--Couplings and Masses #!TO BE FINE_TUNED
# λ=1
# α=1
# β=1
# η=1
# m_ρ=1
# m_ϕ=sqrt(λ*η^2)
# m_ψ=1

# #--Data Recorder
# nsnaps=200

end #module
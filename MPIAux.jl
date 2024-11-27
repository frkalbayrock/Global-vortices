module MPIAux

include("Parameters.jl")
using .Parameters
using MPI
using DelimitedFiles

#!MPI PART
const comm = MPI.COMM_WORLD
const myrank = MPI.Comm_rank(comm)
const nprocs = MPI.Comm_size(comm)
const comm_cart = MPI.Cart_create(comm, nprocs_perdim; periodic=periods, reorder=false)
# const coords_cart = MPI.Cart_coords(comm_cart, myrank)

const buffRecv_ϕ_x = im*zeros(padd,Ny_loc)
const buffRecv_ψ_x =    zeros(padd,Ny_loc)
const buffRecv_Z_x = im*zeros(padd, Nx, Ny_loc, Ny) #this is along x-direction
const buffRecv_ϕ_y = im*zeros(Nx_loc,padd)
const buffRecv_ψ_y =    zeros(Nx_loc,padd)
const buffRecv_Z_y = im*zeros(Nx_loc, Nx, padd, Ny) #this is along y-direction


export chunker
function chunker(rank)

    coords_cart_rank = MPI.Cart_coords(comm_cart, rank)

    lx_loc = Int((Nx_loc * coords_cart_rank[1] - (Nx/2-1)))
    rx_loc = Int(lx_loc + Nx_loc - 1)
    ly_loc = Int((Ny_loc * coords_cart_rank[2] - (Ny/2-1)))
    ry_loc = Int(ly_loc + Ny_loc - 1)

    return lx_loc, rx_loc, ly_loc, ry_loc
end
 



# export  find_neighbours
function find_neighbours()
    #Left-right
    direction=0
    disp=1 
    n_left, n_right = MPI.Cart_shift(comm_cart,direction,disp)
    #Bottom-top
    direction=1
    disp=1
    n_bottom, n_top = MPI.Cart_shift(comm_cart,direction,disp)

return n_left, n_right, n_bottom, n_top
end

const n_left, n_right, n_bottom, n_top = find_neighbours()


#Prepares the local chunk fields to derivative operations
#by updating the paddings built into locak chunks.
export update_Paddings!
function update_Paddings!(ϕ,ψ,Z)
    transfer_Paddings_ϕ_ψ!(ϕ,buffRecv_ϕ_x,buffRecv_ϕ_y)
    transfer_Paddings_ϕ_ψ!(ψ,buffRecv_ψ_x,buffRecv_ψ_y)
    transfer_Paddings_Z!(Z,buffRecv_Z_x,buffRecv_Z_y)
    MPI.Barrier(comm)
end




function transfer_Paddings_ϕ_ψ!(f,buffRecv_x,buffRecv_y)

    #Left to right
    @views MPI.Sendrecv!(f[Nx_loc+padd:end-padd, padd+1:Ny_loc+padd], n_right, 1,  #!!!!!
                    buffRecv_x, n_left, 1, comm)
    @views f[1:padd,padd+1:Ny_loc+padd] .= buffRecv_x


    #Right to left
    @views MPI.Sendrecv!(f[padd+1:2*padd, padd+1:Ny_loc+padd], n_left, 2, #!!!!!!!!
                    buffRecv_x, n_right, 2, comm)
    @views f[Nx_loc+padd+1:end,padd+1:Ny_loc+padd] .= buffRecv_x


    #Bottom to top
    @views MPI.Sendrecv!(f[padd+1:Nx_loc+padd,Ny_loc+padd:end-padd], n_top, 3, #!!!!!!!
                    buffRecv_y, n_bottom, 3, comm)
    @views f[padd+1:Nx_loc+padd, 1:padd] .= buffRecv_y
    

    #Top to bottom
    @views MPI.Sendrecv!(f[padd+1:Nx_loc+padd, padd+1:2*padd], n_bottom, 4,     #!!!!!
                    buffRecv_y , n_top, 4, comm)
    @views f[padd+1:Nx_loc+padd, Ny_loc+padd+1:end] .= buffRecv_y

end






function transfer_Paddings_Z!(Z,buffRecv_x,buffRecv_y)#,n_left, n_right, n_bottom, n_top)

    #Left to right
    @views MPI.Sendrecv!(Z[Nx_loc+padd:end-padd,:, padd+1:Ny_loc+padd,:], n_right, 1,   #!!!!!
                    buffRecv_x, n_left, 1, comm)
    @views Z[1:padd,:, padd+1:Ny_loc+padd,:] .= buffRecv_x


    #Right to left
    @views MPI.Sendrecv!(Z[padd+1:2*padd,:, padd+1:Ny_loc+padd,:], n_left, 2,           #!!!!!
                    buffRecv_x, n_right, 2, comm)
    @views Z[Nx_loc+padd+1:end,:, padd+1:Ny_loc+padd,:] .= buffRecv_x


    #Bottom to top
    @views MPI.Sendrecv!(Z[padd+1:Nx_loc+padd, :, Ny_loc+padd:end-padd,:], n_top, 3,    #!!!!!
                    buffRecv_y, n_bottom, 3, comm)
    @views Z[padd+1:Nx_loc+padd,:,  1:padd,:] .= buffRecv_y
    

    #Top to bottom
    @views MPI.Sendrecv!(Z[padd+1:Nx_loc+padd,:, padd+1:2*padd,:], n_bottom, 4,         #!!!!!
                    buffRecv_y , n_top, 4, comm)
    @views Z[padd+1:Nx_loc+padd,:, Ny_loc+padd+1:end,:] .= buffRecv_y

end

end #module
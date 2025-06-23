module Auxiliary_Routines 

using Distributed
include("Parameters.jl")
using .Parameters
include("IndexMap.jl")
using .IndexMap
using OffsetArrays


#This routine applies the periodic boundary conditions on 2-D lattice
#It takes the lattice point (j,k) and determines nnl (nearest neighbour left) and nnr (nearest neighbour right) on each direction.
#If the point (j,k) is at the boundary the nnl_x,y and nnr_x,y are set accordingly.
#Then it returns nnl_x,y and nnr_x,y.
export pbc2D
function pbc2D(j,k)

    #x-direction
    nnl_x = j-1
    nnr_x = j+1
    nnl_x < lx ?  nnl_x=rx : nothing
    nnr_x > rx ?  nnr_x=lx : nothing

    #y-direction
    nnl_y = k-1
    nnr_y = k+1
    nnl_y < ly ?  nnl_y=ry : nothing
    nnr_y > ry ?  nnr_y=ly : nothing

    return nnl_x, nnr_x, nnl_y, nnr_y 
end


#This routine find the nearest neighbours in x and y direction of a 2D lattice with pbc that is flattened. 
#It is assumed that the flattening is done using row-major order.
export pbc1D
function pbc1D(J)

    #-In x-direction:
    #   nearest neighbour on the right will be  +N away and
    #   nearest neighbour on the left will be -N away
    #   all under mod(N^2)
    if mod(J-N,N2)!=0      #To make sure N^2 is mapped to N^2 not 0
        nnl_x = mod(J-N,N2)    
    else
        nnl_x = N2
    end

    if mod(J+N,N2) !=0     #To make sure N^2 is mapped to N^2 not 0
        nnr_x = mod(J+N,N2)
    else
        nnr_x = N2
    end


    #-In y-direction:
    #   nearest neighbour on the right will be  +1 away and
    #   nearest neighbour on the left will be -1 away
    #   unless they are the last or first element of their block respectively
    #   For those we wrap around inside the block: 
    #   first to last (for the left neighbour) and last to first (for the right neighbour) is mapped.
    if mod(J,N) == 1
        nnl_y = J+N-1
    else
        nnl_y = J-1
    end

    if  mod(J,N) == 0
        nnr_y = J-N+1
    else
        nnr_y = J+1
    end

    return nnl_x, nnr_x, nnl_y, nnr_y
end




#Returns the coordinates of the global ends of the lattice chunk given proc id
export chunker
function chunker(proc_id::Int)
    
    #Cartesian coordinate of the chunk
    i_procs, j_procs = chunk_cart(proc_id)

    #Global coordinates of the ends of the chunks
    #(assumes OffsetArrays)
    lx_p = Int((Nx_loc * i_procs - (Nx/2-1)))
    rx_p = Int(lx_p + Nx_loc - 1)
    ly_p = Int((Ny_loc * j_procs - (Ny/2-1)))
    ry_p = Int(ly_p + Ny_loc - 1)

return lx_p, rx_p, ly_p, ry_p
end




#Returns the Cartesian coordinate of the chunk of a given proc id
function chunk_cart(proc_id::Int)
    i_procs = Int(mod((proc_id-2),nprocs_perdim[1]))
    j_procs = Int((proc_id-2 - i_procs)/nprocs_perdim[1])

return i_procs, j_procs
end




#Periodic boundary conditions on the chunks
#Gets the Cartesian coord of the chunk &
#returns neighbouring coordinates along x and y directions
function chunk_pbc(i::Int,j::Int)

    #x-direction
    nnl = i-1
    nnr = i+1
    nnl < 0 ? nnl=nprocs_perdim[1]-1 : nnl
    nnr > nprocs_perdim[1]-1 ? nnr=0 : nnr
    
    #y-direction
    nnb = j-1
    nnt = j+1
    nnb < 0 ? nnb=nprocs_perdim[2]-1 : nnb
    nnt > nprocs_perdim[2]-1 ? nnt=0 : nnt

return nnl, nnr, nnb, nnt
end




#Gets the Cartesian coordinate of the chunk and gives back proc_id
function chunk_id(i_procs::Int,j_procs::Int)
    proc_id = j_procs * nprocs_perdim[1] + i_procs + 2
return proc_id
end




export find_neighbours
function find_neighbours(proc_id::Int)
    
    #Get Cartesian coordinate of the chunk first
    i_procs, j_procs = chunk_cart(proc_id)
    nnl, nnr, nnb, nnt = chunk_pbc(i_procs,j_procs) #nearest neighbours l:left, r:right, b:bottom, t:top

    #Neighbours
    n_left =   chunk_id(nnl, j_procs)
    n_right =  chunk_id(nnr, j_procs)
    n_bottom = chunk_id(i_procs, nnb)
    n_top =    chunk_id(i_procs, nnt)

return n_left, n_right, n_bottom, n_top
end




#This part is used initizalie the 'RemoteChannel's for data transfers.
#This method is very similar to MPI but using only the Distributed.jl.
#'RemoteChannel's are the references/handles of 'Channel's that can reside on any worker used to store data.
#The following initizaliation is conducted at two steps;
    #1-Each worker creates the receiving data channel(s) locally
    #2-The 'RemoteChannel's/handles are exhange between neigbors.
        #This step requires the master to coordinate the operation since
        #@fetchfrom'ing neigbors' channels with the same name causes small problems. 
        #(since master has no channels residing in it, this doesn't cause a problem over there)
        #Thus, master fetches the handles and distributes to the workers who are neighbors.
#We use a single receiving channel on each worker, where the identification of messages are done with symbols: :from_left etc.
#All the neighbors dump data into this channel with proper tags (the symbols).
#Then receive function sorts things out later to update the proper sites. That way we don't deal with many channels and remotecalls.
export initialize_channels
function initialize_channels()
    
    # Step 1: Each worker creates its own receive buffers locally
    @everywhere workers() const recv_ψ = RemoteChannel(() -> Channel{Tuple{Symbol, Array{Float64,2}}}(4))
    @everywhere workers() const recv_ϕ = RemoteChannel(() -> Channel{Tuple{Symbol, Array{ComplexF64,2}}}(4))
    @everywhere workers() const recv_Z = RemoteChannel(() -> Channel{Tuple{Symbol, Array{ComplexF64,4}}}(4))
    
    # Step 2: Exchange handles (push buffer handles to neighbors who will send into them)
    @sync for pid in workers()
        left, right, bottom, top = find_neighbours(pid)

        @async begin
            #ψ
            ch_ψ = @fetchfrom pid Main.recv_ψ
            #ϕ
            ch_ϕ = @fetchfrom pid Main.recv_ϕ
            #Z
            ch_Z = @fetchfrom pid Main.recv_Z

            # Send these channels to left/right/bottom/top neighbors
            remotecall_wait(left) do
                global Main.send_ψ_right = ch_ψ
                global Main.send_ϕ_right = ch_ϕ
                global Main.send_Z_right = ch_Z
            end
            remotecall_wait(right) do
                global Main.send_ψ_left = ch_ψ
                global Main.send_ϕ_left = ch_ϕ
                global Main.send_Z_left = ch_Z
            end
            remotecall_wait(bottom) do
                global Main.send_ψ_top = ch_ψ
                global Main.send_ϕ_top = ch_ϕ
                global Main.send_Z_top = ch_Z
            end
            remotecall_wait(top) do
                global Main.send_ψ_bottom = ch_ψ
                global Main.send_ϕ_bottom = ch_ϕ
                global Main.send_Z_bottom = ch_Z
            end
        end
    end

end




#Data exchange initiator to fill paddings with updated data from neighboring blocks
export update_Paddings!
function update_Paddings!()

    #Send in to the buffers
    @everywhere  workers() begin
        Main.send_f!()
    end 

    #Collect from the buffers
    @everywhere  workers() begin
        Main.recv_f!()
    end

return nothing
end




#Send Functions - send the borders to paddings of the neighbors.
@everywhere function send_f!()
    #ψ
    put!(send_ψ_left,   (:from_right,  ψ[padd+1:padding_size,  padd+1:Ny_loc+padd]))
    put!(send_ψ_right,  (:from_left,   ψ[Nx_loc+1:end-padd,    padd+1:Ny_loc+padd]))
    put!(send_ψ_bottom, (:from_top,    ψ[padd+1:Nx_loc+padd,   padd+1:padding_size]))
    put!(send_ψ_top,    (:from_bottom, ψ[padd+1:Nx_loc+padd,   Ny_loc+1:end-padd]))
    #ϕ
    put!(send_ϕ_left,   (:from_right,  ϕ[padd+1:padding_size,  padd+1:Ny_loc+padd]))
    put!(send_ϕ_right,  (:from_left,   ϕ[Nx_loc+1:end-padd,    padd+1:Ny_loc+padd]))
    put!(send_ϕ_bottom, (:from_top,    ϕ[padd+1:Nx_loc+padd,   padd+1:padding_size]))
    put!(send_ϕ_top,    (:from_bottom, ϕ[padd+1:Nx_loc+padd,   Ny_loc+1:end-padd]))
    #Z
    put!(send_Z_left,   (:from_right,  Z[padd+1:padding_size, :,  padd+1:Ny_loc+padd, :]))
    put!(send_Z_right,  (:from_left,   Z[Nx_loc+1:end-padd,   :,  padd+1:Ny_loc+padd, :]))
    put!(send_Z_bottom, (:from_top,    Z[padd+1:Nx_loc+padd,  :,  padd+1:padding_size,:]))
    put!(send_Z_top,    (:from_bottom, Z[padd+1:Nx_loc+padd,  :,  Ny_loc+1:end-padd,  :]))
end




#Receive Functions - receive borders of the neigbors to the paddings.
@everywhere function recv_f!()
    # ψ
    for _ in 1:4
        dir, data = take!(recv_ψ)
        if dir == :from_left
            ψ[1:padd, padd+1:Ny_loc+padd] = data
        elseif dir == :from_right
            ψ[Nx_loc+padd+1:end, padd+1:Ny_loc+padd] = data
        elseif dir == :from_bottom
            ψ[padd+1:Nx_loc+padd, 1:padd] = data
        elseif dir == :from_top
            ψ[padd+1:Nx_loc+padd, Ny_loc+padd+1:end] = data
        end
    end

    # ϕ
    for _ in 1:4
        dir, data = take!(recv_ϕ)
        if dir == :from_left
            ϕ[1:padd, padd+1:Ny_loc+padd] = data
        elseif dir == :from_right
            ϕ[Nx_loc+padd+1:end, padd+1:Ny_loc+padd] = data
        elseif dir == :from_bottom
            ϕ[padd+1:Nx_loc+padd, 1:padd] = data
        elseif dir == :from_top
            ϕ[padd+1:Nx_loc+padd, Ny_loc+padd+1:end] = data
        end
    end

    # Z
    for _ in 1:4
        dir, data = take!(recv_Z)
        if dir == :from_left
            Z[1:padd, :, padd+1:Ny_loc+padd, :] = data
        elseif dir == :from_right
            Z[Nx_loc+padd+1:end, :, padd+1:Ny_loc+padd, :] = data
        elseif dir == :from_bottom
            Z[padd+1:Nx_loc+padd, :, 1:padd, :] = data
        elseif dir == :from_top
            Z[padd+1:Nx_loc+padd, :, Ny_loc+padd+1:end, :] = data
        end
    end
end





end #module
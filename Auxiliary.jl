module Auxiliary_Routines 

using Distributed
include("Parameters.jl")
using .Parameters
include("IndexMap.jl")
using .IndexMap
using OffsetArrays
using DistributedArrays
using DistributedArrays.SPMD


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
    if mod(J-N,N^2)!=0      #To make sure N^2 is mapped to N^2 not 0
        nnl_x = mod(J-N,N^2)    
    else
        nnl_x = N^2
    end

    if mod(J+N,N^2) !=0     #To make sure N^2 is mapped to N^2 not 0
        nnr_x = mod(J+N,N^2)
    else
        nnr_x = N^2
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






#Returns the ends of each chunk given the proc_id --(no padding)--
export distChunker
function distChunker(proc_id::Int)
    
    #Cartesian coordinate of the chunk
    i_procs, j_procs = chunk_cart(proc_id)

    #Global coordinates of the ends of the chunks 
    lx_loc = Int((Nx_loc * i_procs +1))
    rx_loc = Int(lx_loc + Nx_loc -1)
    ly_loc = Int((Ny_loc * j_procs +1))
    ry_loc = Int(ly_loc + Ny_loc - 1)

return lx_loc, rx_loc, ly_loc, ry_loc
end



#Returns the ends of each chunk for fields --with "padding"-- given the proc_id
export distChunkerPadd
function distChunkerPadd(proc_id::Int)
    
    #Cartesian coordinate of the chunk
    i_procs, j_procs = chunk_cart(proc_id)

    #Get the sizes of each block
    Nx_chunk = Int(Nx_padd/nprocs_perdim[1])
    Ny_chunk = Int(Ny_padd/nprocs_perdim[2])

    #Global coordinates of the ends of the chunks 
    lx_loc = Int((Nx_chunk * i_procs +1))
    rx_loc = Int(lx_loc + Nx_chunk -1)
    ly_loc = Int((Ny_chunk * j_procs +1))
    ry_loc = Int(ly_loc + Ny_chunk - 1)

return lx_loc, rx_loc, ly_loc, ry_loc
end





#---Transfer Functions
#Data exchange initiator to fill paddings with updated data from neighboring blocks
export update_Paddings!
function update_Paddings!(ϕ,ψ,Z)
    
    update_Padds!(ϕ)
    update_Padds!(ψ)
    update_Padds!(Z)

    # @time begin
    #     #!v3 for updating paddings using DistArray directly
    #     @sync @distributed for p in workers()
    #         update_Padds2!(ψ)
    #         update_Padds2!(ϕ)
    #         update_Padds2!(Z)
    #     end
    # end

end

#!---------TURN ONLY ONE OF THEM FOR TESTING FOR THE BIGGER LATTICE!-----------!#
#!version1 (so far this seems to be faster and more memory efficient)
function update_Padds!(f::DArray{<:Number, 2})

    #Find neighbors
    n_left, n_right, n_bottom, n_top = find_neighbours(myid())

    #!we'll add @views at some point.
    #Left padding
    f[:L][1:padd, padd+1:Ny_loc+padd]            .= @fetchfrom n_left   f[:L][Nx_loc+padd:end-padd, padd+1:Ny_loc+padd]
    #Right padding
    f[:L][Nx_loc+padd+1:end, padd+1:Ny_loc+padd] .= @fetchfrom n_right  f[:L][padd+1:padding_size, padd+1:Ny_loc+padd]
    #Bottom padding
    f[:L][padd+1:Nx_loc+padd, 1:padd]            .= @fetchfrom n_bottom f[:L][padd+1:Nx_loc+padd, Ny_loc+padd:end-padd]
    #Top padding
    f[:L][padd+1:Nx_loc+padd, Ny_loc+padd+1:end] .= @fetchfrom n_top    f[:L][padd+1:Nx_loc+padd, padd+1:padding_size]

return nothing
end


#Updating 
function update_Padds!(Z::DArray{ComplexF64, 4})

    #Find neighbors
    n_left, n_right, n_bottom, n_top = find_neighbours(myid())

    #!we'll add @views at some point.
    #Left padding
    Z[:L][1:padd, :, padd+1:Ny_loc+padd, :]            .= @fetchfrom n_left   Z[:L][Nx_loc+padd:end-padd, :, padd+1:Ny_loc+padd, :]
    #Right padding
    Z[:L][Nx_loc+padd+1:end, :, padd+1:Ny_loc+padd, :] .= @fetchfrom n_right  Z[:L][padd+1:padding_size, :, padd+1:Ny_loc+padd, :]
    #Bottom padding
    Z[:L][padd+1:Nx_loc+padd, :, 1:padd, :]            .= @fetchfrom n_bottom Z[:L][padd+1:Nx_loc+padd, :, Ny_loc+padd:end-padd, :]
    #Top padding
    Z[:L][padd+1:Nx_loc+padd, :, Ny_loc+padd+1:end, :] .= @fetchfrom n_top    Z[:L][padd+1:Nx_loc+padd, :, padd+1:padding_size, :]

return nothing
end



#!version2
# function update_Padds!(f::DArray{<:Number, 2})

#     #Find neighbors
#     n_left, n_right, n_bottom, n_top = find_neighbours(myid())

#     #--Send
#     sendto(n_right,  f[:L][Nx_loc+padd:end-padd, padd+1:Ny_loc+padd]; tag="lr") #Left to right
#     sendto(n_left,   f[:L][padd+1:padding_size, padd+1:Ny_loc+padd];  tag="rl") #Right to left
#     sendto(n_top,    f[:L][padd+1:Nx_loc+padd, Ny_loc+padd:end-padd]; tag="bt") #Bottom to top
#     sendto(n_bottom, f[:L][padd+1:Nx_loc+padd, padd+1:padding_size];  tag="tb") #Top to bottom
    
#     #--Receive
#     f[:L][1:padd, padd+1:Ny_loc+padd]            .= recvfrom(n_left;   tag="lr")
#     f[:L][Nx_loc+padd+1:end, padd+1:Ny_loc+padd] .= recvfrom(n_right;  tag="rl")
#     f[:L][padd+1:Nx_loc+padd, 1:padd]            .= recvfrom(n_bottom; tag="bt")
#     f[:L][padd+1:Nx_loc+padd, Ny_loc+padd+1:end] .= recvfrom(n_top;    tag="tb")


#     barrier(;pids=workers())
# end

# function update_Padds!(Z::DArray{ComplexF64, 4})

#     #Find neighbors
#     n_left, n_right, n_bottom, n_top = find_neighbours(myid())

#     #--Send
#     sendto(n_right,  Z[:L][Nx_loc+padd:end-padd, :, padd+1:Ny_loc+padd, :]; tag="lr") #Left to right
#     sendto(n_left,   Z[:L][padd+1:padding_size, :, padd+1:Ny_loc+padd, :];  tag="rl") #Right to left
#     sendto(n_top,    Z[:L][padd+1:Nx_loc+padd, :, Ny_loc+padd:end-padd,: ]; tag="bt") #Bottom to top
#     sendto(n_bottom, Z[:L][padd+1:Nx_loc+padd, :, padd+1:padding_size, :];  tag="tb") #Top to bottom
    
#     #--Receive
#     Z[:L][1:padd, :, padd+1:Ny_loc+padd, :]            .= recvfrom(n_left;   tag="lr")
#     Z[:L][Nx_loc+padd+1:end, :, padd+1:Ny_loc+padd, :] .= recvfrom(n_right;  tag="rl")
#     Z[:L][padd+1:Nx_loc+padd, :, 1:padd, :]            .= recvfrom(n_bottom; tag="bt")
#     Z[:L][padd+1:Nx_loc+padd, :, Ny_loc+padd+1:end, :] .= recvfrom(n_top;    tag="tb")

#     barrier(;pids=workers())
# end
#!---------TURN ONLY ONE OF THEM FOR TESTING FOR THE BIGGER LATTICE!-----------!#








#!#############################--------------------------Version - 3 -- looks way slower
    # export update_Padds2!
    # function update_Padds2!(ψ::DArray{<:Number, 2})

    #     lx_p, rx_p, ly_p, ry_p = distChunkerPadd(myid())
    #     #Find neighbors
    #     n_left, n_right, n_bottom, n_top = find_neighbours(myid())
    #     # lx_left, rx_left, ly_left, ry_left = distChunkerPadd(n_left)

    #     #!how to do pbc??
    #     nnl = lx_p-padd-1
    #     nnr = rx_p+padd+1
    #     nnb = ly_p-padd-1
    #     nnt = ry_p+padd+1

    #     nnl < 1       && (nnl = (Nx_padd)-padd)
    #     nnr > Nx_padd && (nnr = 1+padd)
    #     nnb < 1       && (nnb = (Ny_padd)-padd)
    #     nnt > Ny_padd && (nnt = 1+padd)


    #     #Left        
    #     # ψ[lx_p:lx_p+padd-1, ly_p+padd:ry_p-padd] .= ψ[nnl:nnl+padd-1, ly_p+padd:ry_p-padd]
    #     localpart(ψ)[1:padd, padd+1:Ny_loc+padd] .= ψ[nnl:nnl+padd-1, ly_p+padd:ry_p-padd]
    #     #Right
    #     # ψ[rx_p-padd+1:rx_p, ly_p+padd:ry_p-padd, t] .= ψ[nnr:nnr+padd-1, ly_p+padd:ry_p-padd, t]
    #     localpart(ψ)[Nx_loc+padd+1:end, padd+1:Ny_loc+padd] .= ψ[nnr:nnr+padd-1, ly_p+padd:ry_p-padd]
    #     #Bottom
    #     # ψ[lx_p+padd:rx_p-padd, ly_p:ly_p+padd-1, t] .= ψ[lx_p+padd:rx_p-padd, nnb:nnb+padd-1, t]
    #     localpart(ψ)[padd+1:Nx_loc+padd, 1:padd] .= ψ[lx_p+padd:rx_p-padd, nnb:nnb+padd-1]
    #     #Top
    #     # ψ[lx_p+padd:rx_p-padd, ry_p-padd+1:rx_p, t] .= ψ[lx_p+padd:rx_p-padd, nnt:nnt+padd-1, t]
    #     localpart(ψ)[padd+1:Nx_loc+padd, Ny_loc+padd+1:end] .= ψ[lx_p+padd:rx_p-padd, nnt:nnt+padd-1]

    # return nothing
    # end

    # export update_Padds2!
    # function update_Padds2!(Z::DArray{ComplexF64, 4})

    #     lx_p, rx_p, ly_p, ry_p = distChunkerPadd(myid())
    #     #Find neighbors
    #     n_left, n_right, n_bottom, n_top = find_neighbours(myid())
    #     # lx_left, rx_left, ly_left, ry_left = distChunkerPadd(n_left)

    #     #!how to do pbc??
    #     nnl = lx_p-padding_size
    #     nnr = rx_p+padd+1
    #     nnb = ly_p-padding_size
    #     nnt = ry_p+padd+1

    #     nnl < 1 && (nnl = (Nx_padd+1)-padding_size)
    #     nnr > Nx_padd && (nnr = 1+padd)
    #     nnb < 1       && (nnb = (Ny_padd+1)-padding_size)
    #     nnt > Ny_padd && (nnt = 1+padd)


    #     #Left        
    #     localpart(Z)[1:padd, :, padd+1:Ny_loc+padd, :] .= Z[nnl:nnl+padd-1, :, ly_p+padd:ry_p-padd, :]
    #     #Right
    #     localpart(Z)[Nx_loc+padd+1:end, :, padd+1:Ny_loc+padd, :] .= Z[nnr:nnr+padd-1, :, ly_p+padd:ry_p-padd, :]
    #     #Bottom
    #     localpart(Z)[padd+1:Nx_loc+padd, :, 1:padd, :] .= Z[lx_p+padd:rx_p-padd, :, nnb:nnb+padd-1, :]
    #     #Top
    #     localpart(Z)[padd+1:Nx_loc+padd, :, Ny_loc+padd+1:end, :] .= Z[lx_p+padd:rx_p-padd, :, nnt:nnt+padd-1, :]

    # return nothing
    # end

#!#############################--------------------------Version - 2








export chunker
function chunker(proc_id::Int)
    
    #Cartesian coordinate of the chunk
    i_procs, j_procs = chunk_cart(proc_id)

    #Global coordinates of the ends of the chunks 
    lx_loc = Int((Nx_loc * i_procs - (Nx/2-1)))
    rx_loc = Int(lx_loc + Nx_loc - 1)
    ly_loc = Int((Ny_loc * j_procs - (Ny/2-1)))
    ry_loc = Int(ly_loc + Ny_loc - 1)

return lx_loc, rx_loc, ly_loc, ry_loc
end





#Gives the Cartesian coordinate of the chunk
function chunk_cart(proc_id::Int)
    i_procs = Int(mod((proc_id-2),nprocs_perdim[1]))
    j_procs = Int((proc_id-2 - i_procs)/nprocs_perdim[1])

return i_procs, j_procs
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





#Gets the Cartesian coordinate of the chunk and gives back proc_id
function chunk_id(i_procs::Int,j_procs::Int)
    proc_id = j_procs * nprocs_perdim[1] + i_procs + 2
return proc_id
end





#Periodic boundary conditions on the chunks
#Gets the Cartesian coord of the chunk &
#returns neighbouring coordinates along x and y directions
#Warning: Chunk coordinate system is a 0-based coordinate system, 
#meaning the the coordinates of the block myid()=2 (which is the first workers()) => (0,0)
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





# #---Transfer Functions
# #   So far this is the only way it works with @fetchfrom that's why it is done not so clever way, on purpose.

# #Data exchange initiator to fill paddings with updated data from neighboring blocks
# # function update_Paddings(t::Int)
# @everywhere workers() function update_Paddings(t::Int)
#         #Find neighbors
#         n_left, n_right, n_bottom, n_top = find_neighbours(myid())
#         #Notation: neighbour: 1=left, 2=right, 3=bottom, 4=top
#         ϕ[1:padd, padd+1:Ny_loc+padd, t]            .= @fetchfrom n_left   getData_ϕ(1,t)
#         ϕ[Nx_loc+padd+1:end, padd+1:Ny_loc+padd, t] .= @fetchfrom n_right  getData_ϕ(2,t)
#         ϕ[padd+1:Nx_loc+padd, 1:padd, t]            .= @fetchfrom n_bottom getData_ϕ(3,t)
#         ϕ[padd+1:Nx_loc+padd, Ny_loc+padd+1:end, t] .= @fetchfrom n_top    getData_ϕ(4,t)

#         ψ[1:padd, padd+1:Ny_loc+padd, t]            .= @fetchfrom n_left   getData_ψ(1,t)
#         ψ[Nx_loc+padd+1:end, padd+1:Ny_loc+padd, t] .= @fetchfrom n_right  getData_ψ(2,t)
#         ψ[padd+1:Nx_loc+padd, 1:padd, t]            .= @fetchfrom n_bottom getData_ψ(3,t)
#         ψ[padd+1:Nx_loc+padd, Ny_loc+padd+1:end, t] .= @fetchfrom n_top    getData_ψ(4,t)

#         Z[1:padd, :, padd+1:Ny_loc+padd, :, t]            .= @fetchfrom n_left   getData_Z(1,t)
#         Z[Nx_loc+padd+1:end, :, padd+1:Ny_loc+padd, :, t] .= @fetchfrom n_right  getData_Z(2,t)
#         Z[padd+1:Nx_loc+padd, :, 1:padd, :, t]            .= @fetchfrom n_bottom getData_Z(3,t)
#         Z[padd+1:Nx_loc+padd, :, Ny_loc+padd+1:end, :, t] .= @fetchfrom n_top    getData_Z(4,t)
# return nothing
# end


# # @everywhere workers() 
# # export getData_ψ
# # function getData_ψ(neighbour::Int,t::Int)
# @everywhere workers() function getData_ψ(neighbour::Int,t::Int) #! for now i have to define them like this for name space issues, 
#                                                                     #! also needed "using Distributed" in this module
#                                                                     #!I wanna make Auxiliary not a module eventually.
#     #Notation: neighbour: 1=left, 2=right, 3=bottom, 4=top
#     if neighbour == 1
#         return @views ψ[Nx_loc+padd:end-padd, padd+1:Ny_loc+padd, t]::SubArray{Float64, 2, Array{Float64, 3}, Tuple{UnitRange{Int64}, UnitRange{Int64}, Int64}, false}
#     elseif neighbour == 2
#         return @views ψ[padd+1:padding_size, padd+1:Ny_loc+padd, t]::SubArray{Float64, 2, Array{Float64, 3}, Tuple{UnitRange{Int64}, UnitRange{Int64}, Int64}, false}
#     elseif neighbour == 3
#         return  @views ψ[padd+1:Nx_loc+padd, Ny_loc+padd:end-padd, t]::SubArray{Float64, 2, Array{Float64, 3}, Tuple{UnitRange{Int64}, UnitRange{Int64}, Int64}, false}
#     elseif neighbour == 4
#         return @views ψ[padd+1:Nx_loc+padd, padd+1:padding_size, t]::SubArray{Float64, 2, Array{Float64, 3}, Tuple{UnitRange{Int64}, UnitRange{Int64}, Int64}, false}
#     else
#         error("Invalid neighbour index")
#     end
# end


# # @everywhere workers() 
# # export getData_ϕ
# # function getData_ϕ(neighbour::Int,t::Int)
# @everywhere workers() function getData_ϕ(neighbour,t)
#     #Notation: neighbour: 1=left, 2=right, 3=bottom, 4=top
#     if neighbour == 1
#         return @views ϕ[Nx_loc+padd:end-padd, padd+1:Ny_loc+padd, t]::SubArray{ComplexF64, 2, Array{ComplexF64, 3}, Tuple{UnitRange{Int64}, UnitRange{Int64}, Int64}, false}
#     elseif neighbour == 2
#         return @views ϕ[padd+1:padding_size, padd+1:Ny_loc+padd, t]::SubArray{ComplexF64, 2, Array{ComplexF64, 3}, Tuple{UnitRange{Int64}, UnitRange{Int64}, Int64}, false}
#     elseif neighbour == 3
#         return @views ϕ[padd+1:Nx_loc+padd, Ny_loc+padd:end-padd, t]::SubArray{ComplexF64, 2, Array{ComplexF64, 3}, Tuple{UnitRange{Int64}, UnitRange{Int64}, Int64}, false}
#     elseif neighbour == 4
#         return @views ϕ[padd+1:Nx_loc+padd, padd+1:padding_size, t]::SubArray{ComplexF64, 2, Array{ComplexF64, 3}, Tuple{UnitRange{Int64}, UnitRange{Int64}, Int64}, false}
#     else
#         error("Invalid neighbour index")
#     end
# end


# # @everywhere workers() 
# # export getData_Z
# # function getData_Z(neighbour::Int,t::Int)
# @everywhere workers() function getData_Z(neighbour,t)
#     #Notation: neighbour: 1=left, 2=right, 3=bottom, 4=top
#     if neighbour == 1
#         return @views Z[Nx_loc+padd:end-padd, :, padd+1:Ny_loc+padd, :, t]::SubArray{ComplexF64, 4, Array{ComplexF64, 5}, Tuple{UnitRange{Int64}, Base.Slice{Base.OneTo{Int64}}, UnitRange{Int64}, Base.Slice{Base.OneTo{Int64}}, Int64}, false}
#     elseif neighbour == 2
#         return @views Z[padd+1:padding_size, :, padd+1:Ny_loc+padd, :, t]::SubArray{ComplexF64, 4, Array{ComplexF64, 5}, Tuple{UnitRange{Int64}, Base.Slice{Base.OneTo{Int64}}, UnitRange{Int64}, Base.Slice{Base.OneTo{Int64}}, Int64}, false}
#     elseif neighbour == 3
#         return @views Z[padd+1:Nx_loc+padd, :, Ny_loc+padd:end-padd, :, t]::SubArray{ComplexF64, 4, Array{ComplexF64, 5}, Tuple{UnitRange{Int64}, Base.Slice{Base.OneTo{Int64}}, UnitRange{Int64}, Base.Slice{Base.OneTo{Int64}}, Int64}, false}
#     elseif neighbour == 4
#         return @views Z[padd+1:Nx_loc+padd, :, padd+1:padding_size, :, t]::SubArray{ComplexF64, 4, Array{ComplexF64, 5}, Tuple{UnitRange{Int64}, Base.Slice{Base.OneTo{Int64}}, UnitRange{Int64}, Base.Slice{Base.OneTo{Int64}}, Int64}, false}
#     else
#         error("Invalid neighbour index")
#     end
# end


#

end #module
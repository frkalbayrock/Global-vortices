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



#---Transfer Functions
#   So far this is the only way it works with @fetchfrom that's why it is done not so clever way, on purpose.

#Data exchange initiator to fill paddings with updated data from neighboring blocks
@everywhere workers() function update_Paddings!(ϕ,ψ,Z)
    
    #We @async data acquisition so we can stack up the processes from each neighbor for time saving.
    #Notation: neighbour: 1=left, 2=right, 3=bottom, 4=top
    @sync begin
        #Left
        @async begin
            left_data   = @fetchfrom n_left     getData_All(1)::@NamedTuple{ϕ::Matrix{ComplexF64}, ψ::Matrix{Float64}, Z::Array{ComplexF64, 4}}
            ϕ[1:padd, padd+1:Ny_loc+padd]       .= left_data.ϕ
            ψ[1:padd, padd+1:Ny_loc+padd]       .= left_data.ψ
            Z[1:padd, :, padd+1:Ny_loc+padd, :] .= left_data.Z
        end

        #Right
        @async begin
            right_data  = @fetchfrom n_right   getData_All(2)::@NamedTuple{ϕ::Matrix{ComplexF64}, ψ::Matrix{Float64}, Z::Array{ComplexF64, 4}}
            ϕ[Nx_loc+padd+1:end, padd+1:Ny_loc+padd]       .= right_data.ϕ
            ψ[Nx_loc+padd+1:end, padd+1:Ny_loc+padd]       .= right_data.ψ
            Z[Nx_loc+padd+1:end, :, padd+1:Ny_loc+padd, :] .= right_data.Z
        end

        #Bottom
        @async begin
            bottom_data = @fetchfrom n_bottom getData_All(3)::@NamedTuple{ϕ::Matrix{ComplexF64}, ψ::Matrix{Float64}, Z::Array{ComplexF64, 4}}
            ϕ[padd+1:Nx_loc+padd, 1:padd]       .= bottom_data.ϕ
            ψ[padd+1:Nx_loc+padd, 1:padd]       .= bottom_data.ψ
            Z[padd+1:Nx_loc+padd, :, 1:padd, :] .= bottom_data.Z
        end

        #Top
        @async begin
            top_Data    = @fetchfrom n_top     getData_All(4)::@NamedTuple{ϕ::Matrix{ComplexF64}, ψ::Matrix{Float64}, Z::Array{ComplexF64, 4}}
            ϕ[padd+1:Nx_loc+padd, Ny_loc+padd+1:end]       .= top_Data.ϕ
            ψ[padd+1:Nx_loc+padd, Ny_loc+padd+1:end]       .= top_Data.ψ
            Z[padd+1:Nx_loc+padd, :, Ny_loc+padd+1:end, :] .= top_Data.Z
        end

    end
return nothing
end



@everywhere workers() function getData_All(neighbour::Int)

    return (
        ϕ = getData_ϕ(neighbour),
        ψ = getData_ψ(neighbour),
        Z = getData_Z(neighbour)
    )

end


# @everywhere workers() 
# export getData_ψ
# function getData_ψ(neighbour::Int,t::Int)
@everywhere workers() function getData_ψ(neighbour::Int) #! for now i have to define them like this for name space issues, 
                                                                    #! also needed "using Distributed" in this module
                                                                    #!I wanna make Auxiliary not a module eventually.
    #Notation: neighbour: 1=left, 2=right, 3=bottom, 4=top
    if neighbour == 1
        return ψ[Nx_loc+1:end-padd,   padd+1:Ny_loc+padd]::Array{Float64, 2}
    elseif neighbour == 2
        return ψ[padd+1:padding_size, padd+1:Ny_loc+padd]::Array{Float64, 2}
    elseif neighbour == 3
        return ψ[padd+1:Nx_loc+padd,  Ny_loc+1:end-padd]::Array{Float64, 2}
    elseif neighbour == 4
        return ψ[padd+1:Nx_loc+padd,  padd+1:padding_size]::Array{Float64, 2}
    else
        error("Invalid neighbour index")
    end
end


# @everywhere workers() 
# export getData_ϕ
# function getData_ϕ(neighbour::Int,t::Int)
@everywhere workers() function getData_ϕ(neighbour::Int)
    #Notation: neighbour: 1=left, 2=right, 3=bottom, 4=top
    if neighbour == 1
        return ϕ[Nx_loc+1:end-padd,   padd+1:Ny_loc+padd]::Array{ComplexF64, 2}
    elseif neighbour == 2
        return ϕ[padd+1:padding_size, padd+1:Ny_loc+padd]::Array{ComplexF64, 2}
    elseif neighbour == 3
        return ϕ[padd+1:Nx_loc+padd,  Ny_loc+1:end-padd]::Array{ComplexF64, 2}
    elseif neighbour == 4
        return ϕ[padd+1:Nx_loc+padd,  padd+1:padding_size]::Array{ComplexF64, 2}
    else
        error("Invalid neighbour index")
    end
end


# @everywhere workers() 
# export getData_Z
# function getData_Z(neighbour::Int,t::Int)
@everywhere workers() function getData_Z(neighbour::Int)
    #Notation: neighbour: 1=left, 2=right, 3=bottom, 4=top
    if neighbour == 1
        return Z[Nx_loc+1:end-padd,    :, padd+1:Ny_loc+padd,  :]::Array{ComplexF64, 4}
    elseif neighbour == 2
        return Z[padd+1:padding_size,  :, padd+1:Ny_loc+padd,  :]::Array{ComplexF64, 4}
    elseif neighbour == 3
        return Z[padd+1:Nx_loc+padd,   :, Ny_loc+1:end-padd,   :]::Array{ComplexF64, 4}
    elseif neighbour == 4
        return Z[padd+1:Nx_loc+padd,   :, padd+1:padding_size, :]::Array{ComplexF64, 4}
    else
        error("Invalid neighbour index")
    end
end


#

end #module
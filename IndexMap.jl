module IndexMap

include("Parameters.jl")
using .Parameters 
using OffsetArrays

export flattenDimension
#The routine to map a 2-dimensional lattice into a 1-dimensional lattice by row-major ordering
#This particular routine maps the "function values" from 2d lattice to 1d lattice.
#"Function values" can be fields, energy density etc.
#It is necesssary for usage of CQC which requires matrices. 
#[CQC in 1d gives a matrix but 2d gives a 4-index object (not a tensor)]
function flattenDimension(f)

    #Check if the incoming field is real or complex for memory purposes(we define f_s accordingly)
    if eltype(f) == Float64
        f_s = zeros(Nx*Ny)          #"s" stands for single
    elseif eltype(f) == ComplexF64
        f_s = im*zeros(Nx*Ny)
    else
        println("Error: given type not float or complex float.")
    end

    for j=lx:rx
        for k=ly:ry
            #Get the flatten ordered index
            J = twoIndexToOne(j,k)
            f_s[J] = f[j,k]
        end
    end

    return f_s
end

export ravelDimension
#This routine maps back the "function values" 1-d (row-major) flattened dimension to the 2-d lattice dimensions.
#This particular routine maps the "function values" from 1d flattened lattice to 2d physical lattice.
#"Function values" can be fields, energy density etc.
function ravelDimension(f)

    #Check if the incoming field is real or complex for memory purposes(we define f_s accordingly)
    if eltype(f) == Float64
        f_t = zeros(Nx,Ny)          #"s" stands for single
    elseif eltype(f) == ComplexF64
        f_t = im*zeros(Nx,Ny)
    else
        println("---> Error: given type not float or complex float.")
        println(typeof(f))
    end
    f_t = OffsetArray(f_t,lx:rx,ly:ry)

    for J=1:N^2
        j,k=oneIndexToTwo(J)
        f_t[j,k] = f[J]
    end

return f_t
end

export mapZTo4Index
#This routine is particular for Z and dZdt which have 4-indices in 2d CQC setting.
#It takes the flattened Z_JK values and maps them to Z_jlkm 
#where 'j' and 'k' are for space coordinates and 'l' and 'm' are the auxillary indices for CQC.
function mapZTo4Index(Z,dZdt)

    Z_t = im*zeros(Nx,Nx,Ny,Ny)
    dZdt_t = im*zeros(Nx,Nx,Ny,Ny)
    Z_t = OffsetArray(Z_t,lx:rx,lx:rx,ly:ry,ly:ry)
    dZdt_t = OffsetArray(dZdt_t,lx:rx,lx:rx,ly:ry,ly:ry)

    for J=1:N^2
        for K=1:N^2
            j,k = oneIndexToTwo(J)
            l,m = oneIndexToTwo(K)
            Z_t[j,l,k,m] = Z[J,K]
            dZdt_t[j,l,k,m] = dZdt[J,K]
        end
    end

return Z_t, dZdt_t
end

#Reducing a 2d lattice into a 1d one using row-major ordering. (Flattenning)
#Takes two index and reduces to one and returns.
function twoIndexToOne(j,k)
    #Bring the indices from lattice coordinates lx:rx to conventional ordering 1:N
    j += Nx/2
    k += Ny/2
    #Map the 2d coordinate to 1d one
    J = Int((j-1)*N + k)
    
    return J
end

#Inverse map from flattened lattice to back to 2d lattice.
#Takes an index from 1d latt and returns the coordinates in the 2d physical lattice.
function oneIndexToTwo(J)
    j = Int(round(J/N,RoundUp) - Nx/2)
    k = Int(mod(J-1,N)+1 - Ny/2)
    #Nx/2 and Ny/2 factors are to shift the indices for the lattice coordinates

return j,k
end



end
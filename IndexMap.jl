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

    # #-Version 1
        # #Check if the incoming field is real or complex for memory purposes(we define f_s accordingly)
        # if eltype(f) == Float64
        #     # f_s = zeros(Nx*Ny)          #"s" stands for single
        #     f_s = Array{Float64,1}(undef, Nx*Ny)    #"s" stands for single
        # elseif eltype(f) == ComplexF64
        #     # f_s = im*zeros(Nx*Ny)
        #     f_s = Array{ComplexF64,1}(undef, Nx*Ny)
        # else
        #     println("Error: given type not float or complex float.")
        # end

        # for j=lx:rx
        #     for k=ly:ry
        #         #Get the flatten ordered index
        #         J = twoIndexToOne(j,k)
        #         f_s[J] = f[j,k]
        #     end
        # end

    #-Version 2 (uses the intrinsic flatten function from Julia)
        f_s = collect(Iterators.flatten(f'))    #! seems to be faster and uses less allocations/memory

return f_s
end




export ravelDimension
#This routine maps back the "function values" 1-d (row-major) flattened dimension to the 2-d lattice dimensions.
#This particular routine maps the "function values" from 1d flattened lattice to 2d physical lattice.
#"Function values" can be fields, energy density etc.
function ravelDimension(f)

    #Check if the incoming field is real or complex for memory purposes(we define f_s accordingly)
    if eltype(f) == Float64
        f_t = OffsetArray(Array{Float64,2}(undef, Nx,Ny),lx:rx,ly:ry)    #"t" stands for "two-index"
    elseif eltype(f) == ComplexF64
        f_t = OffsetArray(Array{ComplexF64,2}(undef, Nx,Ny),lx:rx,ly:ry)
    else
        println("---> Error: given type not float or complex float.")
        println(typeof(f))
    end

    for J=1:N2
        j,k=oneIndexToTwo(J)
        f_t[j,k] = f[J]
    end

return f_t
end




export mapZTo4Index
#This routine is particular for Z and dZdt which have 4-indices in 2d CQC setting.
#It takes the flattened Z_JK (or dZdt_JK) values and maps them to Z_jlkm (or dZdt_jlkm)
#where 'j' and 'k' are for space coordinates and 'l' and 'm' are the auxillary indices for CQC.
function mapZTo4Index(ZordZ)

    if eltype(ZordZ) == Float64
            ZordZ_t = Array{Float64,4}(undef, Nx,Nx,Ny,Ny)    #"s" stands for single
        elseif eltype(ZordZ) == ComplexF64
            ZordZ_t = Array{ComplexF64,4}(undef, Nx,Nx,Ny,Ny)
        else
            println("Error: given type not float or complex float.")
    end
    # ZordZ_t = im*zeros(Nx,Nx,Ny,Ny) #t stands for "two-index" for each given index (i.e. j,k <- J)#!
    # ZordZ_t = OffsetArray(ZordZ_t,lx:rx,lx:rx,ly:ry,ly:ry)#! ***** turn this on when you remove the bottom red warning.

    for J=1:N2
        for K=1:N2
            j,k = oneIndexToTwo(J)
            l,m = oneIndexToTwo(K)
            #! REMOVE THIS ONCE WE ARE DONE USING 4-INDEX OMEGA'S also make sure you turn on **** above
            j= Int(j+Nx/2)
            k= Int(k+Ny/2)
            l= Int(l+Nx/2)
            m= Int(m+Ny/2)
            #! 
            ZordZ_t[j,l,k,m] = ZordZ[J,K]
        end
    end

return ZordZ_t
end




export mapZTo2Index
function mapZTo2Index(ZordZ)

    #Version 1 (direct mapping)
        #     ZordZ_f = Array{ComplexF64,2}(undef, Nx*Ny,Nx*Ny)  #f stands for "four"
        #     for j=lx:rx
        #         for k=ly:ry
        #             for l=lx:rx
        #                 for m=ly:ry
        #                     J = twoIndexToOne(j,k)
        #                     K = twoIndexToOne(l,m)
        #                     ZordZ_f[J,K] = ZordZ[j,l,k,m]
        #                 end
        #             end
        #         end
        #     end
        # return ZordZ_f

    #Version 2 (uses the flattenDimension function)
        #Map the space coordinates
        Z_spaceFlatten = Array{ComplexF64,3}(undef, Nx*Ny,Nx,Ny)
        Z_spaceFlatten = OffsetArray(Z_spaceFlatten, 1:Nx*Ny,lx:rx, ly:ry)
        for l=lx:rx
            for m=ly:ry
                Z_spaceFlatten[:,l,m] .= flattenDimension(ZordZ[:,l,:,m])
            end
        end

        #Map the auxillary coordinates
        Z_Flatten = Array{ComplexF64,2}(undef, Nx*Ny,Nx*Ny)
        for J=1:Nx*Ny
            Z_Flatten[J,:] .= flattenDimension(Z_spaceFlatten[J,:,:])
        end

return Z_Flatten
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
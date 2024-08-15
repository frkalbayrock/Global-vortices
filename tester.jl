# include("Auxiliary.jl")
# using .Auxiliary_Routines 
# using OffsetArrays


#Testing pbc1D
# N=3
# for J=1:N^2
#     nnl_x, nnr_x, nnl_y, nnr_y = pbc1D(J,N)
#     println("---J = ",J,"---")
#     println("nnr_x: ",nnr_x," // nnl_x: ", nnl_x, " // nnr_y: ", nnr_y, " // nnl_y: ",nnl_y)
# end


#Testing shifting indices from j=1,N to j=lx,rx
# N=4
# Nx=N
# Ny=N
# for J=1:N^2
#     j = Int(round(J/N,RoundUp)) - Nx/2
#     k = Int(mod(J,N)) - Ny/2
#     println("j: ",j, " // k: ", k)
# end 



#Multiple line caluclation tester
# x=1
#     +1

#     println(x)


# #Testing local/global variables for for loop
# ϕ = zeros(-1:1) .+1
# println(ϕ)

# # for




#Test broadcasting when defining an OffsetArray



# #Test OffsetArray for very large arrays. #!don't work need more ram!
# N=100000
# Nx=Ny=N
# lx=Int.(-Nx/2+1)
# rx=Int.(Nx/2)
# ly=lx
# ry=rx

# @time begin
#     Z_t = im*zeros(Nx,Ny,2)
#     dZdt_t = im*zeros(Nx,Ny,2)
# end

# @time begin
#     Z_t = im*zeros(Nx,Ny,2)
#     dZdt_t = im*zeros(Nx,Ny,2)
#     Z_t = OffsetArray(Z_t,lx:rx,ly:ry,0:1)
#     dZdt_t = OffsetArray(dZdt_t,lx:rx,ly:ry,0:1)
# end



#Testing multi-threading
using Base.Threads

i = Threads.Atomic{Int}(0);
ids = zeros(4);
old_is = zeros(4);
Threads.@threads for id in 1:4
    old_is[id] = Threads.atomic_add!(i, id)
    ids[id] = id
end
println(old_is)
println(i[])
println(i)
println(ids)
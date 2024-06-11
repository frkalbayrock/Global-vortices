include("auxiliary.jl")
using .Auxiliary_routines


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

x=1
    +1

    println(x)
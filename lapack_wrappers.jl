#-Tridiagonalize ---> dsytrd from LAPACK
"""
    Julia doesn't have a warpper for the tridiagonalizetion function dsytrd from LAPACK.
    Here we call it directly by using ccal and define our own wrapper.
    The output is a SymTridiagonal matrix.
"""
const BlasInt = LinearAlgebra.LAPACK.BlasInt
"""
    tridiagonalize_with_Q(A::Matrix{Float64})

Compute the symmetric tridiagonal reduction of `A` using LAPACK's `dsytrd!`,
and return both the tridiagonal matrix `H` and the orthogonal matrix `Q`
such that `A ≈ Q * H * Q'`.
"""
function sym_tridiagonalize!(A::Matrix{Float64})
    n = size(A, 1)
    d = zeros(Float64, n)
    e = zeros(Float64, n - 1)
    tau = zeros(Float64, n - 1)
    info = Ref{BlasInt}(0)

    # === Workspace query for DSYTRD ===
    lwork = BlasInt(-1)
    work = zeros(Float64, 1)
    ccall((LAPACK.@blasfunc(dsytrd_), LAPACK.liblapack), Cvoid,
        (Ref{UInt8}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt},
         Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64},
         Ref{BlasInt}, Ref{BlasInt}),
        'U', BlasInt(n), A, BlasInt(n), d, e, tau, work, lwork, info)

    lwork = BlasInt(work[1])
    work = zeros(Float64, lwork)

    # === Actual DSYTRD call ===
    ccall((LAPACK.@blasfunc(dsytrd_), LAPACK.liblapack), Cvoid,
        (Ref{UInt8}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt},
         Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64},
         Ref{BlasInt}, Ref{BlasInt}),
        'U', BlasInt(n), A, BlasInt(n), d, e, tau, work, lwork, info)
    info[] == 0 || error("dsytrd! failed with info=$(info[])")

    # === Now generate Q from reflectors ===
    lwork = BlasInt(-1)
    work = zeros(Float64, 1)
    ccall((LAPACK.@blasfunc(dorgtr_), LAPACK.liblapack), Cvoid,
        (Ref{UInt8}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt},
         Ptr{Float64}, Ptr{Float64}, Ref{BlasInt}, Ref{BlasInt}),
        'U', BlasInt(n), A, BlasInt(n), tau, work, lwork, info)

    lwork = BlasInt(work[1])
    work = zeros(Float64, lwork)

    ccall((LAPACK.@blasfunc(dorgtr_), LAPACK.liblapack), Cvoid,
        (Ref{UInt8}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt},
         Ptr{Float64}, Ptr{Float64}, Ref{BlasInt}, Ref{BlasInt}),
        'U', BlasInt(n), A, BlasInt(n), tau, work, lwork, info)
    info[] == 0 || error("dorgtr! failed with info=$(info[])")

    Q = copy(A)
    # H = Tridiagonal(e, d, copy(e))
    H = SymTridiagonal(d, e)
return H, Q
end




#-Tridiagonal Eigen Solver ---> DTEDC(Divide and conquer method)
"""
    stedc_manual(T::SymTridiagonal{Float64}; compute_vectors::Bool=true)

Compute eigenvalues (and optionally eigenvectors) of a symmetric tridiagonal
matrix T using LAPACK's DSTEDC (divide-and-conquer).

Returns (w, Z) where `w` is a Vector{Float64} of length n,
and `Z` is n×n (or empty if compute_vectors=false).
"""
function stedc_manual(T::SymTridiagonal{Float64}; compute_vectors::Bool=true)
    n = length(T.dv)
    if n == 0
        return Float64[], Array{Float64}(undef, 0, 0)
    end

    # LAPACK expects arrays that it may overwrite
    d = copy(T.dv)
    # for DSTEDC E should be length n-1
    e = copy(T.ev)

    compz = compute_vectors ? 'I' : 'N'   # 'I' -> compute eigenvectors of tridiagonal
    ldz = max(1, n)
    Z = compute_vectors ? Matrix{Float64}(undef, ldz, n) : Matrix{Float64}(undef, 1, 1)

    info = Ref{BlasInt}(0)

    # workspace query: lwork = -1, liwork = -1
    work = Vector{Float64}(undef, 1)
    iwork = Vector{BlasInt}(undef, 1)

    ccall((LAPACK.@blasfunc(dstedc_), LAPACK.liblapack), Cvoid,
        (Ref{UInt8}, Ref{BlasInt}, Ptr{Float64}, Ptr{Float64},
         Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt},
         Ptr{BlasInt}, Ref{BlasInt}, Ref{BlasInt}),
        compz, BlasInt(n), d, e,
        Z, BlasInt(ldz), work, BlasInt(-1),
        iwork, BlasInt(-1), info)

    if info[] != 0
        error("Workspace-query dstedc returned INFO = $(info[]).")
    end

    # read suggested workspace sizes
    lwork = max(1, BlasInt(Int(work[1])))
    liwork = max(1, iwork[1])

    work = Vector{Float64}(undef, lwork)
    iwork = Vector{BlasInt}(undef, liwork)

    # actual call
    ccall((LAPACK.@blasfunc(dstedc_), LAPACK.liblapack), Cvoid,
        (Ref{UInt8}, Ref{BlasInt}, Ptr{Float64}, Ptr{Float64},
         Ptr{Float64}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt},
         Ptr{BlasInt}, Ref{BlasInt}, Ref{BlasInt}),
        compz, BlasInt(n), d, e,
        Z, BlasInt(ldz), work, lwork,
        iwork, liwork, info)

    if info[] != 0
        error("LAPACK dstedc failed with INFO = $(info[]).")
    end

    # d now contains eigenvalues in ascending order, Z contains eigenvectors (if requested)
    w = copy(d)
return w, compute_vectors ? Z[:, 1:n] : Array{Float64}(undef, 0, 0)
end

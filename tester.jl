# include("Auxiliary.jl")
# using .Auxiliary_Routines 
# using OffsetArrays


#!Testing pbc1D
# N=3
# for J=1:N^2
#     nnl_x, nnr_x, nnl_y, nnr_y = pbc1D(J,N)
#     println("---J = ",J,"---")
#     println("nnr_x: ",nnr_x," // nnl_x: ", nnl_x, " // nnr_y: ", nnr_y, " // nnl_y: ",nnl_y)
# end


#!Testing shifting indices from j=1,N to j=lx,rx
# N=4
# Nx=N
# Ny=N
# for J=1:N^2
#     j = Int(round(J/N,RoundUp)) - Nx/2
#     k = Int(mod(J,N)) - Ny/2
#     println("j: ",j, " // k: ", k)
# end 



#!Multiple line caluclation tester
# x=1
#     +1

#     println(x)


# #!Testing local/global variables for for loop
# ϕ = zeros(-1:1) .+1
# println(ϕ)

# # for




#!Test broadcasting when defining an OffsetArray



# #!Test OffsetArray for very large arrays. #!don't work need more ram!
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



# #!Testing multi-threading
# using Base.Threads

# i = Threads.Atomic{Int}(0);
# ids = zeros(4);
# old_is = zeros(4);
# Threads.@threads for id in 1:4
#     old_is[id] = Threads.atomic_add!(i, id)
#     ids[id] = id
# end
# println(old_is)
# println(i[])
# println(i)
# println(ids)

# #!Testing synchronization
# using Distributed
# addprocs(5)

# @everywhere function my_function()
#     println("Worker $(myid()) is doing some work")
#     sleep(rand())  # Simulate work
#     if myid()==2
#         A=rand(1000000000)
#     end
#     println("Worker $(myid()) finished work")
# end

# @everywhere function my_function2()
#     println("Worker $(myid()) is doing some work-2")
#     sleep(rand())  # Simulate work
#     println("Worker $(myid()) finished work-2")
# end

# # @sync begin
# #     for p in workers()
# #         @async remotecall_wait(my_function, p)
# #     end
# # end

# # @sync begin 
# #     for p in workers()
# #         @async remotecall_wait(my_function2, p)
# #     end
# # end


# @everywhere  workers() my_function()

# # @everywhere workers() my_function2()


# println("All workers have finished executing my_function")







#!Another synchronization test
# using Distributed
# addprocs(4)

# #Setup
# @everywhere begin
#     const mychannel = Channel{Int}(1)
    
#     function send_stuff()
#         println("Process ", myid(), " sending...")
#         put!(mychannel, 42)
#         println("Process ", myid(), " finished sending")
#     end

#     function recv_stuff()
#         println("Process ", myid(), " receiving...")
#         x = take!(mychannel)
#         println("Process ", myid(), " received ", x)
#     end
# end




# # ##---1
# # # First all processes will try sending:
# # @everywhere workers() send_stuff()

# # # Now all processes will try receiving:
# # @everywhere workers() recv_stuff()


    

# ##----2
# # Phase 1: synchronized sending
# @sync for p in workers()
#     @async remotecall_wait(send_stuff, p)
# end

# # Phase 2: synchronized receiving
# @sync for p in workers()
#     @async remotecall_wait(recv_stuff, p)
# end



# #!Understanding RemoteChannels
# #!1-
# using Distributed
# addprocs(2)   # Add 2 worker processes

# @everywhere using Distributed

# # Create a RemoteChannel that lives on process 2:
# ch = RemoteChannel(()->Channel{Int}(10), 2)

# # Now process 1 can do:
# put!(ch, 42)

# # And process 2 (or any process) can do:
# x = take!(ch)   # returns 42


# #!2-
# using Distributed
# addprocs(4)

# @everywhere using Distributed

# # Create a work queue on process 1
# workqueue = RemoteChannel(()->Channel{Int}(10+4), 1)

# # Fill it with work
# for i in 1:10
#     put!(workqueue, i)
# end

# # Add stop signals
# for i in 1:nworkers()
#     put!(workqueue, 109)
# end

# # Workers process work
# @everywhere function worker_main(ch)
#     while true
#         job = take!(ch)
#         if job == 109
#             println("Worker $(myid()) exiting.")
#             break
#         else
#             println("Worker $(myid()) processing job $job")
#         end
#     end
# end

# # Launch workers
# @sync for pid in workers()
#     @async remotecall_wait(worker_main, pid, workqueue)
# end



#!One-Dimensional Point-to-Point Exchange (2 neighbors)
#!v1- notworking properly
# using Distributed

# addprocs(4)

# @everywhere using Distributed

# # Set up neighbors in 1D ring
# const pids = workers()
# const neighbors = Dict(pid => (
#     left = pids[mod1(findfirst(==(pid), pids) - 1, length(pids))],
#     right = pids[mod1(findfirst(==(pid), pids) + 1, length(pids))]
# ) for pid in pids)

# # Each process creates its own receive channels
# @everywhere const recv_ψ = Dict{Symbol, RemoteChannel}()

# for pid in pids
#     left_ch =  RemoteChannel(()->Channel{Vector{Float64}}(1), pid)
#     right_ch = RemoteChannel(()->Channel{Vector{Float64}}(1), pid)

#     # Send these to neighbors
#     remotecall_wait(neighbors[pid][:right]) do
#         global recv_ψ
#         recv_ψ[:left] = right_ch
#     end
#     remotecall_wait(neighbors[pid][:left]) do
#         global recv_ψ
#         recv_ψ[:right] = left_ch
#     end
# end

# # Define ψ and halo send/recv
# @everywhere workers() begin
#         const ψ = fill(myid(), 10)
# end

# @everywhere begin
#     function send_ψ!()
#         put!(recv_ψ[:left],  ψ[2:3])   # send left halo
#         put!(recv_ψ[:right], ψ[8:9])   # send right halo
#     end

#     function recv_ψ!()
#         ψ[1:2]       = take!(recv_ψ[:left])
#         ψ[end-1:end] = take!(recv_ψ[:right])
#     end
# end


# # Check result
# println("BEFORE: ")
# for p in pids
#     # println("PID $p → ", fetch(@spawnat p ψ))
#     anan = @fetchfrom p ψ
#     println("PID $p → ", anan)
    
# end

# # Run exchange
# @sync for p in pids
#     @async remotecall_wait(send_ψ!, p)
# end

# @sync for p in pids
#     @async remotecall_wait(recv_ψ!, p)
# end

# # Check result
# println("AFTER: ")
# for p in pids
#     println("PID $p → ", fetch(@spawnat p ψ))
# end



# #!v2- WORKS!
# using Distributed
# addprocs(4)

# @everywhere using Distributed

# # Worker list and 1D periodic neighbors
# const pids = workers()
# const neighbors = Dict(pid => (
#     left = pids[mod1(findfirst(==(pid), pids) - 1, length(pids))],
#     right = pids[mod1(findfirst(==(pid), pids) + 1, length(pids))]
# ) for pid in pids)



# # Create RemoteChannels only on the receiver side
# const recv_ψ_buffers = Dict{Tuple{Int,Int}, RemoteChannel}()

# for dest_pid in pids
#     for (dir, src_pid) in pairs(neighbors[dest_pid])
#         recv_ψ_buffers[(src_pid, dest_pid)] = RemoteChannel(() -> Channel{Vector{Float64}}(1), dest_pid)
#     end
# end



# # Send correct handles to each worker
# for pid in pids
#     left  = neighbors[pid][:left]
#     right = neighbors[pid][:right]

#     ch_left  = recv_ψ_buffers[(left, pid)]
#     ch_right = recv_ψ_buffers[(right, pid)]

#     remotecall_wait(pid) do
#         global recv_ψ = Dict(
#                     :from_left  => ch_left,
#                     :from_right => ch_right
#                 )
#     end
# end

# # Each worker defines its ψ
# @everywhere workers() begin
#     const ψ = fill(myid(), 10)
# end

# # Define send/receive functions
# @everywhere workers() function recv_ψ!()
#     ψ[1:2]       = take!(recv_ψ[:from_left])
#     ψ[end-1:end] = take!(recv_ψ[:from_right])
# end

# # Define `send_ψ!` on the master (global dictionary)
# @everywhere  workers() function send_ψ!(ch_left, ch_right)
#     put!(ch_left,  ψ[2:3])   # send to left neighbor (they receive from right)
#     put!(ch_right, ψ[8:9])   # send to right neighbor (they receive from left)
# end

# # Show initial state
# println("\nBEFORE:")
# for p in pids
#     println("PID $p → ", @fetchfrom p ψ)
# end

# # Launch sends (use correct channels from global table)
# @sync for pid in pids
#     left  = neighbors[pid][:left]
#     right = neighbors[pid][:right]
#     ch_L  = recv_ψ_buffers[(pid, left)]
#     ch_R  = recv_ψ_buffers[(pid, right)]
#     @async remotecall_wait(pid) do
#         send_ψ!(ch_L, ch_R)
#     end
# end

# # Launch receives
# @sync for pid in pids
#     @async remotecall_wait(() -> recv_ψ!(), pid)
# end

# # Show final state
# println("\nAFTER:")
# for p in pids
#     println("PID $p → ", @fetchfrom p ψ)
# end



# #!v3- trying to simplify --failed....
# using Distributed
# addprocs(4)

# @everywhere using Distributed

# # Worker list and 1D periodic neighbors
# const pids = workers()
# @everywhere const pids = $pids
# @everywhere const neighbors = Dict(pid => (
#     left = pids[mod1(findfirst(==(pid), pids) - 1, length(pids))],
#     right = pids[mod1(findfirst(==(pid), pids) + 1, length(pids))]
# ) for pid in pids)



# # Create RemoteChannels only on the receiver side
# const recv_ψ_buffers = Dict{Tuple{Int,Int}, RemoteChannel}()

# for dest_pid in pids
#     for (dir, src_pid) in pairs(neighbors[dest_pid])
#         recv_ψ_buffers[(src_pid, dest_pid)] = RemoteChannel(() -> Channel{Vector{Float64}}(1), dest_pid)
#     end
# end

# # Send correct handles to each worker
# for pid in pids
#     left  = neighbors[pid][:left]
#     right = neighbors[pid][:right]

#     ch_left  = recv_ψ_buffers[(left, pid)]
#     ch_right = recv_ψ_buffers[(right, pid)]

#     remotecall_wait(pid) do
#     global recv_ψ = Dict(
#                 :from_left  => ch_left,
#                 :from_right => ch_right
#             )
#     end
# end


# @everywhere workers() recv_ψ_buffers = Dict{Tuple{Int,Int}, RemoteChannel}()
# # @everywhere workers() send_ψ_buffers = Dict{Tuple{Int,Int}, RemoteChannel}()

# @everywhere workers() begin

#     n_left  = neighbors[myid()][:left]
#     n_right = neighbors[myid()][:right]

#     # if myid()==3
#     #     @show n_left
#     #     @show n_right
#     # end

#         recv_ψ_buffers[(n_left,  myid())]  = RemoteChannel(
#                                                 ()->Channel{Array{Float64,1}}(1), myid()
#                                             )
#         recv_ψ_buffers[(n_right, myid())] = RemoteChannel(
#                                                 ()->Channel{Array{Float64,1}}(1), myid()
#                                             )
#         recv_ψ_buffers[(myid(), n_left)]  = RemoteChannel(
#                                                 ()->Channel{Array{Float64,1}}(1), n_left
#                                             )
#         recv_ψ_buffers[(myid(), n_right)] = RemoteChannel(
#                                                 ()->Channel{Array{Float64,1}}(1), n_right
#                                             )
    
    
#     # if myid()==2 || myid()==3
#     #     @show recv_ψ_buffers[(n_left,  myid())]
#     #     @show recv_ψ_buffers[(n_right, myid())]
#     #     @show recv_ψ_buffers[(myid(), n_left)]
#     #     @show recv_ψ_buffers[(myid(), n_right)]
#     # end
    

#     # send_ψ_buffers[(myid(), n_left)]  = RemoteChannel(
#     #                                                 ()->Channel{Array{Float64,2}}(1), myid()
#     #                                             )
#     # send_ψ_buffers[(myid(), n_right)] = RemoteChannel(
#     #                                                 ()->Channel{Array{Float64,2}}(1), myid()
#     #                                             )
# end


# @everywhere begin 
#     if myid()==2 || myid()==3
#         @show recv_ψ_buffers
#     end
# end

# # Each worker defines its ψ
# @everywhere workers() begin
#     const ψ = fill(myid(), 10)
# end

# # # Define send/receive functions
# # @everywhere workers() function recv_ψ!()
# #     ψ[1:2]       = take!(recv_ψ[:from_left]) #
# #     ψ[end-1:end] = take!(recv_ψ[:from_right])
# # end

# #!
# # Define send/receive functions
# @everywhere workers() function recv_ψ!()
#     ψ[1:2]       = take!(recv_ψ_buffers[(n_left, myid())]) #from_left
#     ψ[end-1:end] = take!(recv_ψ_buffers[(n_right, myid())]) # from_right
# end



# # Define `send_ψ!` on the master (global dictionary)
# @everywhere  workers() function send_ψ!()
#     n_left  = neighbors[myid()][:left]
#     n_right = neighbors[myid()][:right]
#     put!(recv_ψ_buffers[(myid(), n_right)], ψ[2:3])   # send to left neighbor (they receive from right)
#     put!(recv_ψ_buffers[(myid(), n_left)],  ψ[8:9])   # send to right neighbor (they receive from left)
# end

# # Show initial state
# println("\nBEFORE:")
# for p in pids
#     println("PID $p → ", @fetchfrom p ψ)
# end

# # # Launch sends (use correct channels from global table)
# # @sync for pid in pids
# #     left  = neighbors[pid][:left]
# #     right = neighbors[pid][:right]
# #     ch_L  = recv_ψ_buffers[(pid, left)]
# #     ch_R  = recv_ψ_buffers[(pid, right)]
# #     @async remotecall_wait(pid) do
# #         send_ψ!(ch_L, ch_R)
# #     end
# # end

# @sync for pid in pids
#     @async remotecall_wait(() -> send_ψ!(), pid)
# end

# println("ANAN!")

# # # Launch receives
# # @sync for pid in pids
# #     @async remotecall_wait(() -> recv_ψ!(), pid)
# # end

# # Show final state
# println("\nAFTER:")
# for p in pids
#     println("PID $p → ", @fetchfrom p ψ)
# end









# #!v4 - Defining everything on the workers.
# using Distributed
# addprocs(4)

# @everywhere using Distributed

# const pids = workers()
# const neighbors = Dict(pid => (
#     left = pids[mod1(findfirst(==(pid), pids) - 1, length(pids))],
#     right = pids[mod1(findfirst(==(pid), pids) + 1, length(pids))]
# ) for pid in pids)

# function run_anan()
#     # Step 1: Each worker creates its own receive buffers locally
#     @sync for pid in pids
#         @async remotecall_wait(pid) do
#             global recv_ψ = Dict{Symbol, RemoteChannel}()
#             recv_ψ[:from_left]  = RemoteChannel(() -> Channel{Vector{Float64}}(1))
#             recv_ψ[:from_right] = RemoteChannel(() -> Channel{Vector{Float64}}(1))
#         end
#     end

#     # Step 2: Exchange handles (push buffer handles to neighbors who will send into them)
#     @sync for pid in pids
#         left  = neighbors[pid][:left]
#         right = neighbors[pid][:right]

#         @async begin
#             ch_from_left  = @fetchfrom pid recv_ψ[:from_left]
#             ch_from_right = @fetchfrom pid recv_ψ[:from_right]

#             # Send these channels to left/right neighbors
#             remotecall_wait(left) do
#                 global send_ψ_right = ch_from_left  # left sends to my left-facing buffer
#             end

#             remotecall_wait(right) do
#                 global send_ψ_left = ch_from_right  # right sends to my right-facing buffer
#             end
#         end
#     end

#     # Step 3: Define data and send/recv functions
#     @everywhere workers() begin
#         const ψ = fill(myid(), 10)

#         function send_ψ!()
#             put!(send_ψ_left,  ψ[2:3])
#             put!(send_ψ_right, ψ[8:9])
#         end

#         function recv_ψ!()
#             ψ[1:2]       = take!(recv_ψ[:from_left])
#             ψ[end-1:end] = take!(recv_ψ[:from_right])
#         end
#     end

#     # Step 4: Run exchange
#     println("\nBEFORE:")
#     for p in pids
#         println("PID $p → ", @fetchfrom p ψ)
#     end

#     @sync for p in pids
#         @async remotecall_wait(() -> send_ψ!(), p)
#     end

#     @sync for p in pids
#         @async remotecall_wait(() -> recv_ψ!(), p)
#     end

#     println("\nAFTER:")
#     for p in pids
#         println("PID $p → ", @fetchfrom p ψ)
#     end
# end

# @time run_anan()



#!
# using Distributed
# addprocs(4)

# @everywhere using Distributed

# const pids = workers()
# @everywhere const pids = $pids
# @everywhere const neighbors = Dict(pid => (
#     left = pids[mod1(findfirst(==(pid), pids) - 1, length(pids))],
#     right = pids[mod1(findfirst(==(pid), pids) + 1, length(pids))]
# ) for pid in pids)


# function run_anan()
#     # Step 1: Each worker creates its own receive buffers locally
#     @everywhere workers() const recv_ψ = Dict{Symbol, RemoteChannel}()
#     @everywhere workers() begin
#         recv_ψ[:from_left]  = RemoteChannel(() -> Channel{Array{Float64}}(1), myid())
#         recv_ψ[:from_right] = RemoteChannel(() -> Channel{Array{Float64}}(1), myid())
#     end



#     # Step 2: Exchange handles (push buffer handles to neighbors who will send into them)
#     @everywhere workers() begin
#         left  = neighbors[myid()][:left]
#         right = neighbors[myid()][:right]

#             # Get the receive buffers from the neighbors
#             const send_ψ_right = @fetchfrom right recv_ψ[:from_left]
#             const send_ψ_left  = @fetchfrom left  recv_ψ[:from_right]
#     end


#     # Step 3: Define data and send/recv functions
#     @everywhere workers() begin
#         const ψ = fill(myid(), 10)

#         function send_ψ!()
#             put!(send_ψ_left,  ψ[2:3])
#             put!(send_ψ_right, ψ[8:9])
#         end

#         function recv_ψ!()
#             ψ[1:2]       = take!(recv_ψ[:from_left])
#             ψ[end-1:end] = take!(recv_ψ[:from_right])
#         end
#     end

#     # Step 4: Run exchange
#     println("\nBEFORE:")
#     for p in pids
#         println("PID $p → ", @fetchfrom p ψ)
#     end

#     @sync for p in pids
#         @async remotecall_wait(() -> send_ψ!(), p)
#     end

#     @sync for p in pids
#         @async remotecall_wait(() -> recv_ψ!(), p)
#     end

#     println("\nAFTER:")
#     for p in pids
#         println("PID $p → ", @fetchfrom p ψ)
#     end
# end

# @time run_anan()




# #!2D Version
# using Distributed
# addprocs(4)

# @everywhere using Distributed

# const padd=2
# const pids = workers()
# const dims = (2, 2)  # Adjust based on number of workers (e.g. 2x2 grid for 4 workers)
# @assert prod(dims) == length(pids) "Number of workers must match grid size"

# # Arrange workers in row-major order
# const pid_grid = reshape(pids, dims)


# # Build neighbor dictionary
# const neighbors = Dict{Int, NamedTuple{(:left, :right, :top, :bottom), NTuple{4, Int}}}()

# for j in 1:dims[2], i in 1:dims[1]
#     pid = pid_grid[i, j]
#     neighbors[pid] = (
#         left   = pid_grid[mod1(i - 1, dims[1]), j],
#         right  = pid_grid[mod1(i + 1, dims[1]), j],
#         top    = pid_grid[i, mod1(j - 1, dims[2])],
#         bottom = pid_grid[i, mod1(j + 1, dims[2])]
#     )
# end



# function run_anan()
#     # Step 1: Each worker creates its own receive buffers locally
#     @sync for pid in pids
#         @async remotecall_wait(pid) do
#             global recv_ψ = Dict{Symbol, RemoteChannel}()
#             recv_ψ[:from_left]   = RemoteChannel(() -> Channel{Array{Float64,2}}(1))
#             recv_ψ[:from_right]  = RemoteChannel(() -> Channel{Array{Float64,2}}(1))
#             recv_ψ[:from_bottom] = RemoteChannel(() -> Channel{Array{Float64,2}}(1))
#             recv_ψ[:from_top]    = RemoteChannel(() -> Channel{Array{Float64,2}}(1))
#         end
#     end

#     # Step 2: Exchange handles (push buffer handles to neighbors who will send into them)
#     @sync for pid in pids
#         left   = neighbors[pid][:left]
#         right  = neighbors[pid][:right]
#         bottom = neighbors[pid][:bottom]
#         top    = neighbors[pid][:top]

#         @async begin
#             ch_from_left  = @fetchfrom pid recv_ψ[:from_left]
#             ch_from_right = @fetchfrom pid recv_ψ[:from_right]
#             ch_from_bottom = @fetchfrom pid recv_ψ[:from_bottom]
#             ch_from_top = @fetchfrom pid recv_ψ[:from_top]

#             # Send these channels to left/right/bottom/top neighbors
#             remotecall_wait(left) do
#                 global send_ψ_right = ch_from_left
#             end
#             remotecall_wait(right) do
#                 global send_ψ_left = ch_from_right
#             end
#             remotecall_wait(bottom) do
#                 global send_ψ_top = ch_from_bottom
#             end
#             remotecall_wait(top) do
#                 global send_ψ_bottom = ch_from_top
#             end
#         end
#     end

#     # Step 3: Define data and send/recv functions
#     @everywhere workers() begin
#         const ψ = fill(myid(), (8,8))

#         function send_ψ!()
#             put!(send_ψ_left,   ψ[3:4,3:6])
#             put!(send_ψ_right,  ψ[5:6,3:6])
#             put!(send_ψ_bottom, ψ[3:6,3:4])
#             put!(send_ψ_top,    ψ[3:6,5:6])
#         end

#         function recv_ψ!()
#             ψ[1:2,3:6]       = take!(recv_ψ[:from_left])
#             ψ[end-1:end,3:6] = take!(recv_ψ[:from_right])
#             ψ[3:6,1:2]       = take!(recv_ψ[:from_bottom])
#             ψ[3:6,end-1:end] = take!(recv_ψ[:from_top])
#         end
#     end

#     # Step 4: Run exchange
#     println("\nBEFORE:")
#     for p in pids
#         println("PID $p → ", @fetchfrom p ψ)
#     end

#     @sync for p in pids
#         @async remotecall_wait(() -> send_ψ!(), p)
#     end

#     @sync for p in pids
#         @async remotecall_wait(() -> recv_ψ!(), p)
#     end

#     println("\nAFTER:")
#     for p in pids
#         println("PID $p → ", @fetchfrom p ψ)
#     end
# end

# @time run_anan()


# #!Async example ona  single process.
# #this improves the wait time significantly!
# #basically it announces all pid's at the same time and don't wait for one to complete one at a time.
# using Distributed
# addprocs(4)

# @everywhere begin
#     function slow_value()
#         sleep(1)  # simulate delay
#         return myid()
#     end
# end

# # Serial fetch
# @time begin
#     results = []
#     for pid in workers()
#         push!(results, @fetchfrom pid slow_value())
#     end
# end

# # Async fetch
# @time begin
#     results = Channel(length(workers()))
#     @sync for pid in workers()
#         @async put!(results, @fetchfrom pid slow_value())
#     end
# end









#!LazyArrays tests
using BenchmarkTools

N = 500
ϕ = rand(ComplexF64, N, N)
ψ = rand(N, N)
Z = rand(N, N)

mρ² = 1.0
α   = 0.5
β   = 0.3

function naive!(out, ϕ, ψ, Z, mρ², α, β)
    out .= (mρ² .+ α .* abs2.(ϕ) .+ β .* ψ.^2) .* Z
end

out = similar(Z)
@btime naive!(out, ϕ, ψ, Z, mρ², α, β);


using LazyArrays

function lazy!(out, ϕ, ψ, Z, mρ², α, β)
    expr = @~ (mρ² .+ α .* abs2.(ϕ) .+ β .* ψ.^2) .* Z
    copyto!(out, expr)   # evaluate lazily into out
end

out2 = similar(Z)
@btime lazy!(out2, ϕ, ψ, Z, mρ², α, β);

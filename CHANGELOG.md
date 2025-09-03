# CHANGELOG


### v1.5.0 - Vectorization for Speed
We found a way to vectorize the leap_forward() function for time evolution. This gave a huge increase in performance especially for larger lattices. The implentation and tests are done on both small and really large (cluster size) lattices. 

All the changes listed:
#### Energy.jl:
- Unncessary packages loaded for this module were removed. (Profile and Pprof) These can be added if needed for testing but they need to be removed when the testing is done since they increase the memory usage significantly.

#### IC.jl:
- Unncessary packages loaded for this module were removed. (InteractiveUtils) These can be added if needed for testing but they need to be removed when the testing is done since they increase the memory usage significantly.

#### Evolution.jl:
- Macro version of leap_forward and the accompanying macros were removed and moved to the bottom (commented out) in the case of a need of future referral.
- Half_step vectorized version is the active version now. The looped version was moved to the bottom (commented out) in the case of a need for a future referral.
- New vectorized version of leap_forward is added along with the accompanying flux functions. j and k defined as the ranges used in these functions (this is needed since we are working with padded arrays --need to refer to the bulk) for easy reading. 
- In flux_Z!() since the ϕ and ψ arrays are defined as 2D arrays and this vectorization requires all the arrays to be the same size, we define the 4D array M which involves reshaping ϕ and ψ calculation into 4D arrays. This now can be used in a vectorized calculation with 4D Z array.
- Old original looped version of leap_forward is moved to the bottom (commented out) in the case of a need of future referral.


#### main.jl:
- Some @allocated check are removed since we are done with testing.
- For the purpose of fully vectorized time evolution we require the fluxes to be arrays too since, they are used in calculations mixed with arrays. For that purpose we define flux and meanSqrRho as arrays that are const and preallocated in the beginning.


Hopefully we won't need a next step as the cluster tests show promising performance. The only way I can think of improving is maybe the newly added vectorized functions which I don't plan on working on them extensively. We now focus on the testing the physics!





### v1.4.0 - Memory Issue Fixed 
The problem of usage of immense amounts of memory was related to accidently sending full sqrt(Ω) arrays into workers and chunking them at the workers. This full array copying also meant memory usage increases with number of workers. This update fixes that. (We might still have some work to do to improve memory usage;  we fixed the main issue but maybe we can find more improvements)

All the changes listed:
#### Evolution.jl:
- Some unnecessary packages are no longer loaded. 
- Macro is fixed for future testing. It was written for a serial code with time coordinate on the fields. That is now fixed.
- Vectorized half_step!() is commented out and the for-loop-version is now online as we observe better memory usage. And both version are updated/simplified for better readability.
- ϕ_flux and ψ_flux are now stored in a namedtuple with names ϕ and ψ and called with flux.ϕ and flux.ψ.
- Z_flux is returned explicitly flux_Z()


#### IC.jl:
- Some unnecessary packages are no longer loaded.
-  We added the expected types for the "@fetchfrom" of global ϕ and ψ. (Looks like it was skipped the last update)
- Chunking Z is overhauled. Before, we were setting Z on workers with interpolating inv_S_Ωzero_f and S_Ωzero_f which caused copying the whole matrices onto the workers. This was the main reason why we were using a lot of memory. Now we make sure to chunk on the master and then send them to the workers. 
- The inversion operation of the Ω is now done manually rather than using the inverse() function.


#### main.jl:
- All the ComplexF64 arrays are now defined with zeros(Complexf64,Nx,Ny) rather than creating with zeros and then muptilying with imaginary "i".

Next step is to check the allocations and memory usage for the multiple evolution steps.


### v1.3.0 - Upgraded Update Padding Function with Asynchronous Calls
The update_Paddings() function is modified to work with @async to improve performance. The goal was not just to use async calls but also reduce the number of @fetchfrom calls to reduce the remote call overhead. We basically use a bigger messages to send but we do it in a single call now. 

All the changes listed:
#### Auxiliary.jl:
- update_Paddings!() is changed to make less @fetchfrom calls. Instead of getting the data from a particular neighbor one field at a time, now we get them all at the same time as "NamedTuple". Each name in this tuple corresponds to a field. We assign to each field using ..._data.ϕ, ..._data.ψ or ..._data.Z. 
- getData_all() function is created to handle getting all the padd data from a neighbor with one call. The collecting padd from different fields are done on the worker side now.
And since we also @async this operation while left neighbor is preparing/collecting the data we don't wait for it to finish since @async makes sure to call the other neighbors as well. 

#### main.jl:
- BLAS thread count is manually set using the number of cores available. This is also printed on the console later with information about the run. Until now, we weren't aware of the autmatic multi-threading of the BLAS in the background. 


### v1.2.0 - Type Inference Reduction and Performance Update // Symmetric Decleration For Matrix Operations
This update is all about improving the perfomance. We've tested some allocations on the workers (and on master) and we have seen it was (at least on compilation) type inference allocations were dominating it. It is unclear how to test the workers without the compilation so since the results show a lot of type inference we have tested all the functions with "@code_warntype" and fixed all the type stability issues. 
We have also changed how the matrix operations work. Since the matrices we have at initial conditions are all symmetric we define them like that for performance improvement. This way Julia can use the faster LAPACK routines that are for symmetric matrices that reduce the number of operations.
On top of these we have few small changes as well (like "N^2 -> N2") as explained in detail below.

All the changes listed: 
#### Auxiliary.jl:
- All the N^2 are replaced with precomputed N2.
- update_Paddings!() -> update_Paddings!(ϕ,ψ,Z) for consistentcy, readability and possible slight improvement even though it is not really needed.
- Whole update_Paddings!() business slightly changed. It used to be done via "@views" and "SubArray" statements for specifying the type for type stability. Now, we changed that to "Array" statements on "getData_f()"s and we added that type specifiers also on the update_Paddings!(). 
- The second part of the one above is actually part of a general addition. To improve type stability, we added the expected types for all the "@fetchfrom"s. Seems to be improving the type stability.

#### Constrainst.jl:
- All the N^2 are replaced with precomputed N2.

#### Energy.jl:
- ZED -> ZED_gl since now we treat that as a new field rather than defining everytime we call the energy_calculation(). This also means we changed energy_calculation() -> energy_calculation!() as now we are updating "ZED" as a field rather than "return"ing it.
- Initialization totalE = 0 -> totalE = 0. for type stability.
- All "@fetchfrom"s got type statements around them for type stability.

#### Evolution.jl:
- Similar updates as above ones: N^2 -> N2, update_Paddings!() -> update_Paddings!(ϕ,ψ,Z), @fetchfrom got type statements.
- Functions that update fields only got "return nothing" at the end for type stability.

#### IC.jl:
- All the N^2 are replaced with precomputed N2.
- Sparsity tools added. (Checking how many zeros we have in Z) We are keeping these for future testing and for adding Sparse Matrices later.
- Biggest change we have here is the improvement of the calculation of sqrt of matrices. After we define "Ω" matrix we specify it to be a "Symmetric" matrix. Instead of using the given sqrt() for this operation we also use the manual diagonalization path to calculate the sqrt of the matrices. This way Julia knows to use the faster LAPACK routines designed for symmetric matrices for performance improvement.
- All @fetchfrom's got type statements around them for type stability.
- For calculating "Ω" we start with flattenning "ϕ" and "ψ". We have removed @views from those statements as we were already sending the full arrays into the function.

#### IndexMap.jl:
- In this module, some functions needed OffsetArrays but the offsettings were done after the arrays defined. Now we define directly as OffsetArrays. (when needed)
- All the N^2 are replaced with precomputed N2.

#### main.jl:
- Leaving some testing allocation tools for later testing. 
- "ZED" is now defined as a new field on the workers rather than defining everytime when it was needed. This should reduce unnecessary allocations.

#### Parameters.jl:
- N2 is defined from N^2 for precalculation, since this is used a lot in the code.

#### tester.jl:
- We have added a lot of test functions for "RemoteChannel" testing and learning. We have few working ones in there for implementation to the code. We'll be adding another version as  another branch with this new method. "RemoteChannel"s method should build a network of buffers like MPI method directly from Julia/Distributed for consistent networking rather than openning and closing a connection between processes everytime "@fetchfrom" is called. To be tested for performance comparison later.


Next update for this branch will be trying to "@async" "@fetchfrom" calls for improvement. It is unknown at this stage if this is even possible. 





### v1.1.2 - Global Variables and Small Fixes Update
The main change is that the function that determines the neighbours of the chunks is now only called once in the beginning and n_i's are now defined as global constants to avoid repeated calls every time update_Paddings() is called. 

All the changes listed: 

#### Auxiliary.jl:
- Added comments to some functions that lacked before.
- Small variable name change in chunker() to be more clear.
- Repositioned few functions in the order of their call flow.
- update_Paddings() no longer calls find_neigbours(); neighbour are now determined as constants in the beginning of the code.

#### Evolution.jl: 
- macros and functions that only run on the workers are now defined only on the workers.
- updateForNextStep(): now removed the fields as arguments as they are not updated anyway in the leapFrog method.
- fluxes_ϕ_ψ() arguments no longer have specific ranges as the full arrays already being sent. (Thus, we can remove the @views too.)

#### IC.jl:
- When chunking Z, the ends of the chunks were being defined as globals before. We just added "local" in the beginning to make them local variables. Better way of doing this is to put that part into a function which we wrote the function below but there are few issues with it. Couldn't figure out a way to do that yet. 
- As mentioned in the previous point, ic_Z!() is now defined but commented out as it is not working properly (yet). 

#### main.jl:
- The local fields and their time derivatives that are defined on the chunks are now initialized as global constants.
- find_neighbours() is now called in the beginning to determine the neighbours for each workers chunk once rather than calling it everytime we call update_Paddings(). Corresponding n_i's are global constants now.

#### Parameters.jl:
- dx^2 and dy^2 are now precalculated as new constants as dx2 and dy2 to performance on tight loops.


In the next update, it is planned to implement the results of @code_warntype analysis. This will eliminate the type inference allocations. (This one is experimentally done on the serial code and seen lots of improvements already. The plan is to test it on the parallel code and make the necessary changes for improvements on the worker allocations.) 



### v1.1.1 - Small Fixes and Beautification
We have found few possible bugs in getData_f() functions. Also some realignment, organization is done to make the code more readable.

All the changes listed: 

##### Auxiliary.jl:
- In update_Paddings() and getData_f() the array indices are aligned for easy reading.
- There were some extra commas for n=1 and n=2 lines of the getData_ψ() in the arguments of "ψ". These are removed.
- For all getData_f() for n=1 and n=3 lines the data requested didn't scale with derivative order. The issue was for n=1 for example, Nx_loc+padd:end-padd is used for the x-coordinate range which was not correct in general. For our case it worked as this code uses only the nearest neighbors derivative calculation but for higher order derivatives this wouldn't work. Nx_loc+padd:end-padd was targeting only a single column to be copied/sent since Nx_loc+padd = end-padd. If we needed more columns like for higher order derivatives this wouldn't cover those. These statements are changed with the correct expression Nx_loc+1:end-padd which runs on the correct range no matter the order of the derivative. Similar is done for the y-coordinate as well at lines n=3 for each  getData_f() function.

##### Energy.jl:
- We finally added the + λη^4/4 term to the potential energy. This was always done at the run level never applied to the commits; now we finally fixed this.

##### IC.jl:
- Just some realignment for easy reading.

##### main.jl:
- We finally added the command line arguments to create the directories where the data is going to be written. Surprisingly, this was always done manually when needed until now because of being extremely lazy.





### v1.1.0
Tested some new methods; @inline and macros. In evolution module, the fluxes are written with macros or @inline'd. 

All the changes listed:

Evolution.jl:
- Added flux macros in the beginning of the module.
- half_step!() function is updated with ".+" operation which was just "+" before. This was tested earlier on the serial code and we proved that this helped with the vectorization.
- Both flux functions are now @inline'd explicitly both on function definition and at the caller. 

main.jl:
- Added just some benchmarking tools and commands.

Testing (on cluster) shows inconclusive results. However, @inline overall seem to be improving. I'll keep the @inline'd version.

### v1.0.1
Found a way to get rid of the time coordinate only from the fields f (not dfdt's). Could improve memory allocations and overall performance.

All the changes listed:

Auxiliary.jl:
- updatePaddings(), getData_f() no longer have "t" as an argument and inside every "t" is removed.
- In getData_f() we also changed the expected field definitions to be returned, there is one less argument on all of them.

Energy.jl:
- update_Paddings() no longer need argument which is deleted here.
- All the fields inside energy_calculation() now have only space coordinates; all the time dependencies are removed.
- Small cosmetic changes and some extra stuff removed.

Evolution.jl:
- update_Paddings() no longer need argument which is deleted here.
- All the explicit time dependencies from the fields are removed in time_evolution(), half_step(), leap_forward(), flux_Z(), updateForNextStep().
- Small cosmetic changes.

IC.jl: 
- All the explicit time dependencies from the fields are removed in initialConditions(), ic_ϕ,ψ(), renormalization(), zeroPointEnergy(). 
- Some cosmetic changes and some extra stuff removed.

main.jl:
- Definitions of the fields ϕ, ψ and Z were changed to remove the time coordinate.
- Writing the informaiton on a file part is updated as other DArray version of the code with println()'s instead of writedlm().
- 

### v1.0.0
Upgrading the version number only. No changes have been made here.

### v0.3.10
Small changes for parallel-Julia version of the code. Mainly FindVortex module is updated for vortex positions recording. 

All the changes listed  :

Energy.jl:  
- energy() is now defined energy!()

Evolution.jl:  
- vortex_finder() is used during evolution. It now records every snap time rather than just at the end. And the positions are recorded right after into the same file; first vortices written on a line and then the anti-vortices written on the next line. So odd numbered lines are for vortices and evens are for anti-vortices.
- All the data files used during evolution is closed.

FindVortex.jl:  
- vortex_finder(): vortex_pos and anti_votex_pos defined to record the positions of vortex/anti-vortex and returns these vectors to the calling program.

main.jl:  
- Initial Z energy density (ZED) is nto recorded into a separate file called "initial_ZED.dat".


### v0.3.6
All the changes listed:  

Auxiliary.jl:
- For the parallel functions' arguments, the types are added for performance improvement. 
- update_Padddings() returns nothing now. Before it was accidently returning the last line which had some performance hit.
- getData functions' if statements covers everything now; else added with error statement.

Evolution.jl:
- 

IC.jl:
- Z_gl and dZdt_gl is removed from arguments and local Z's are now calcualted directly from the omega matrices.
- Fluctuations are added for the field ϕ
- Small synhtatic updates.
- S_Ωzero, inv_S_Ωzero are now created first with zeros earlier, then updated with sqrt()'s etc.

main.jl:
- Z_gl and dZdt_gl is no longer created. We set local Z's directly on the workers now using 4-index omega matrices.
- Vortex finder is called after time_evolution()

Parameters.jl:
- lx, rx are defined with Int() rather than Int.16() for performance reasons.
- dt_half is defined here for performance.

Newly added modules
- FindVortex.jl : finds the vortices based on winding number calculated at each plaquette.

Next stage: 
It seems we still have performance issues; possibly due to the way the parallelization is done. Cluster runs show that there is an enormous read and write operations happening. The testing shows that it is due to the parallelization macros. Needs further looking into. 
DistributedArrays might be the only option way forward.

### v0.3.5
Various small improvements and changes on the parallel-julia code. 

All the changes listed:
Auxiliary.jl:
- 2*padd is replaced with padding_size for improved performance. (Even slighest improvements for parallel seems to improve performance)
Energy.jl:
- ZED_gl is now (see updates for main.jl below) feeded to the function as an argument instead defining every time energy() is called and returned.
Evolution.jl:
- ZED_gl is now (see updates for main.jl below) feeded to the function as an argument instead defining every time energy() is called and returned.
- Some confusion in the indices of local variables is sorted out. We now write them explicitly as padd+1:padd+Nx_loc
- Small changes in the running indices of dϕdt, dψdt and dZdt (see main.jl changes below)  
IC.jl:
- Small syntax changes
- zeroPointEnergy() is updated to define Z_partial ath the fetch call not before.
main.jl:
- Global field arrays are define as const global julia variables outside of the run_ev() function. We'll continue testing if this gives significant improvements or not.
- Paddings are removed from the time derivative fields dϕdt,dψdt and dZdt for memory relief. dZdt with the paddings gives a significant memory allocation on paper; about N^3. As a result we also updated all the routines with the time derivatives to run for loops with correct indices. For loops haven't been changed just the dummy indicis in time derivatives are shifted by the amount of padding; i.e. j-> j-padd.
- Time coordinates from the global fields are removed. There is no need since we evolve the local ones. 
- ZED_gl is defined as another field and sent as an argument to the function that were recording or returning it. That way we no longer defined ZED at each recording iteration from scratch. 



### v0.3.2
Z_gl is removed from the timeEvolution() which was put there only for constraints and conserved quantities check.

### v0.3.1
Parallelization is completed. The energy graphs for parallel vs serial code matches exactly.
Only the -4-index Z/2-index ϕ,ψ- notation is used throughout for ease of use of parallelization. 
This might be updated to use the -2-index Z/1-index ϕ,ψ- notation for performance purposes later. 

All the changes listed:
Auxiliary.jl:
- Energy routine is removed from this file. 
- For parallelization all the global to local mappings and transfer of data routines are added.
    - chunker(): finding the (end) coordinates of arrays of local chunks in the global lattice based on their proc id.
    - chunk_cart(): gives the cartesian coordinate of a local chunk based on proc id in the domain decomposition.
    - find_neighbours(): finding the neighbouring chunks in the global domain.
    - chunk_id(): gives the proc id based on the cartesian coordinate of the chunk in the domain; opposite of chunk_cart().
    - update_Paddings(): update the paddings of the local arrays that share boundary information with neighbours. 
    - getData_ϕ,_ψ,_Z(): collects the requested data from a given neighbour.
Evolution.jl:
- Move a time-step part is updated to be used on the workers. update_Paddings() added after the first half-step!() for derivatives in the next step.
- Synchronization is optimized using the master-worker topology. Need to sync everything before update_Paddings() for correct padding updates 
    and also update_Paddings() is synched across workers thus the derivatives are taken once the transfer is over.
- Take a snap part now has a step to collect the field information from the workers (local ϕ,ψ) to master (global ϕ,ψ) for recording.
IC.jl:
- Complete overhaul: ϕ,ψ are now first set in another function sent to all workers thus they are locally defined to begin with.
    ϕ,ψ are then collected to their global field arrays for Ω matrix calculations.
- Global Z's are defined using CQC initial conditions as normal for now. Then they are sent to local chunks.
    - This could be avoided if one can write the initial conditions locally directly. (will be updated in a later version)
        Even with this update Ω matrices will have to be defined on the master globally for matrix operations -but at least removing Z_gl will improve memory overall.
- Renormalization and zero-point energy calculations are now done at the worker that has the corner point (xmax,ymax).
IndexMap.jl:
- mapZto4Index() is under renovation. This is updated to accomodate 4-index Ω matrices This cannot correctly map (offsets and everything) Z to 4-index anymore.
Parameters.jl
- Parameters necessary for parallelization are added.
main.jl
- Global and local arrays are defined on the corresponding procs. 
- For now everything uses 4-index Z/2-index ϕ,ψ notation.
- Definitions of functions have been changed. We no longer need to send in the arrays (mostly).
- Whole main body is defined in a function called run_ev() to avoid defining global (in the julia memory sense) variables.

Newly added modules
- Energy.jl:
    - This is moved from the Auxiliary.jl to here with an additional function.
    - Total energy is now calculated at local chunks and summed over later when needed on the master.
    - ZED is also locally determined and combined into a global ZED when requested.

Next stage:
- Memory and ease of use optimizations
- Finding defects with winding number calculations
- Possibly re-parallelize with the single index notation for better perfomance (need some speed testing first)

### v0.2
Serial code is completed up to finding strings routine.
The code takes given inital conditions, find energy at any stage, time evolve the system by solving the differential equations of the theory and take a snapshot at given intervals. 

All the changes listed:
Auxiliary.jl:
- @views is used wherever possible
- gEϕψ and gEZ are now calculated forward and backward and then average is taken
- Zero-point energy calculation no longer uses the corner of the lattice at every time step:
    Now zPE calculated in the beginning (using the corner) and passed into energy().
IC.jl
- Full physical initial conditions are set for all the fields
- @views is use whenever possible.
- zeroPointEnergy() created for energy calculation
- renormalization() updated with intrinsic Julia function
IndexMap.jl
- flattenDimension() updated with intrinsic function
- ravelDimension() updated for better performance; temporary array is created without allocation
- mapZTo4Index() updated to handle both Z and dZdt (still not memory optimized**)
- mapZto2Index() function created for reserve mapping of Z
main.jl
- New modules are included
- New benchmark and profiling tools are added
- Z and dZdt are implicitly defined for better memory
- "Information about the run" part added
    This includes writing info to a text file as well as to the console
- Initial conditions are now being recorded
Parameters.jl
- All the parameters are defined as const global variables which gave an enormouse memory optimization along with performance.

Newly added modules
- Constrains.jl: Checks the contraint and conserved quantities of CQC for stability purposes (like checking constant energy)
- Evolution.jl: Time evolution of the system by numerically solving the differential equations given by the theory
                This is done in both leap-frog method and Crank-Nicolson method

###### COMPLETE DEBUGGING OF THE ENTIRE CODE IS DONE DUE TO NOT CONSTANT ENERGY. SEARCH FOR THE REASON PUSHED US FOR COMPLETE TESTING OF THE FULL CODE. THE REASON AT THE END IS FOUND TO BE WRONG PHYSICAL INTUITION FOR ZERO-POINT ENERGY.


Next stage:
- Parallelization (CPU)
- Finding defects with winding number calculations




### v0.1
Initial conditions and energy routines are completed
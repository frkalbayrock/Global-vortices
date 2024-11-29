# CHANGELOG

#### v4.0.1 - DArray mixed with Julia Parallel Update
Trying a novel method where we mix the DArray implementation with the julia parallelization macros.
Still work in progress. So far Z dependencies turned off just to make the proof of concept work. Just initial conditions for ϕ and ψ are set and evolved for a sinle time step. So far the concept works. 

All the changes listed:

Auxiliary.jl:
- Various commen-outs nad distChunker() renamed to be chunker().

Evolution.jl:
- Other than a one step evolution everything turned off. 
- All the functions redefined inside and for their arguments to avoid Z and dZdt.
- The rest of the changes for the evolution are the changes to make to concept work. 

IC.jl:
- Similar changes mentioned in the Evolution.jl.

main.jl:
- Similar changes to all above. Nothing extra.



#### v0.3.10 - DistributedArrays Update
DistributedArrays are implemented fully using "@sync @distributed" concept of parallelization.

All the changes listed:

Auxiliary.jl:
- distChunkerPadd() created to handle the lattice with paddings -- now that we are using DistArrays the field arrays need to come with paddings on the full lattice because chunks are not defined individually; there is only one field array. (Apart from f_gl arrays which are for now kept for recording purposes)
- updatePaddings!() is defined. It calls the individual update_Padds!() functions.
- update_Padds!() is defined. There are two definitions one where we send in a 3-indexed DArray (for ϕ and ψ) and the second one is 5-indexed DArray for Z. This function is called inside a @sync @distributed concept. This way we can send in the full DArray "f" rather than "localpart(f)" since the caller is the master proc rather than workers(). This function then takes the DArray as an argument and using the "localpart(f)" and "@fetchfrom" updates the padding of the chunks of the DArray.
- In addition, this whole updating paddings operation is done in another way as a version 2. This version tries to use the DArray directly (at least for the half of the operations) using the global indices for each chunk. What it tries to do is that updating localparts of field arrays by calling the DArray directly with the necessary global indices for that chunk. 
It is abandoned as it looks like it is quite slower than the version 1 described above. 
- All the routines for parallel-Julia version are commented out.

Energy.jl:
- Definitions of energy() ang energy_calculation() has changed to send in the DistArrays.
- update_Paddings() function is called for updating the paddings and total energy is calculated using @sync @distributed with reduction operation (+) for summation of totalE_loc's.
- ZED_loc definition inside energy_calculation() is removed since we define ZED now as a new field which we define them as a DistArray.

Evolution.jl:
- time_evolve() and all other functions' that use fields definitions are updated to receive DArrays.
- Snapshot interval and openning files parts are updated slightly; some comments and some adjustments for println statements.
- Recording of the fields part is updated slighty updated for DArray construction. We have a new function to determine the global indices of local chunks. Another change is the fecthing from localpart()'s.
- Vortex checker is added in snapping part.
- All evolution functions are written with @sync @distributed construction and localpart() function for each field. 

FindVortex.jl:
- Vortex finder has a commented out piece to use directly the DArray rather than ϕ_gl which is defined for recording purposes which might be changed later if we can find a way to just use the padded DArrays rather than f_gl's.

IC.jl:
- Cleaned out the unnecessary junk.
- renormalization() and zeroPointEnergy() are updated to work with DArrays. Instead of @fetchfrom'ing we use DArray and just call what we want like a serial code.

main.jl:
- f_gl's are defined without offset. They are still defined for recording purposes for the time being. We could just record the DArrays but the fields have paddings and recording paddings seems useless. We can test the speed of this vs using f_gl's later. 
- ZED is defined as a DArray as well but doesn't have a global version since it doesn't have a padding.
- Recording the information about lattice/parameters part is updated; writedlm() removed, using just println().
- initial_ZED is recorded.

Next:
We need to test speed recording f_gl() version vs DArray recording (either by removing the paddings before recording or removing paddings in the plotting)
Also is there maybe a way to do all of these with DArray without using the @sync @distributed construct? Like the one they do in DArray documentation but optimized for memory.



#### v0.3.9 - DistributedArrays
Complete new version of parallelization using DistributedArrays.jl package
---Still work in progress---

All the changes listed:

Auxiliary.jl:
- distChunker() created. Finds the (end) coordinates of arrays of local chunks in the global lattice based on their proc id. Same as chunker() but without the "Offset".

FindVortex.jl:
- Updated for recording the positions of vortices and anti-vortices and returns them.

IC.jl:
- localpart(f) is used to let the initial conditions for ϕ and ψ. ϕ_gl and ψ_gl is defined for matrix calculations of initial Z's. Later Z is initialized chunk-wise with the "global" omega matrices.
- Overall indices is changed since DistibutedArrays doesn't support Offset'ting.

Main.jl:
- All the arrays are now defined with indices starting 1. Offset's are now turned off.

Parameters.jl:
- New padded lattice sizes defined for DistributedArray's: Nx_padd,Ny_padd.


#### v.0.3.8
MPI Code updated --> Z_gl and dZdt_gl are no longer defined. Should improve memory usage.

#### v0.3.6
All the changes listed

Auxiliary.jl:
- For the parallel functions' arguments, the types are added for performance improvement. 
- update_Padddings() returns nothing now. Before it was accidently returning the last line which had some performance hit.
- getData functions' if statements covers everything now; else added with error statement.

Evolution.jl:

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

#### v0.3.5
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



#### v0.3.2
Z_gl is removed from the timeEvolution() which was put there only for constraints and conserved quantities check.

#### v0.3.1
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

#### v0.2
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




#### v0.1
Initial conditions and energy routines are completed
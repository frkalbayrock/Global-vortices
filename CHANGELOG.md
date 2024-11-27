# CHANGELOG

#### v2.0.2
MPI version gets the FindVortex module. The vortex_finder() is now called every time step. 

#### v2.0.1
Found a way to get rid of the time coordinate only from the fields f (not dfdt's). Could improve memory allocations and overall performance.

All the changes listed:

Energy.jl:
- update_Paddings() no longer need argument which is deleted here.
- All the fields inside energy_calculation() now have only space coordinates; all the time dependencies are removed.

Evolution.jl:
- update_Paddings() no longer need argument which is deleted here.
- All the explicit time dependencies from the fields are removed in time_evolution(), half_step(), leap_forward(), flux_Z(), updateForNextStep().

IC.jl:
- All the explicit time dependencies from the fields are removed in initialConditions(), ic_ϕ,ψ(), renormalization(), zeroPointEnergy().
- Some cosmetic changes and some extra stuff removed.

main.jl:
- Definitions of the fields ϕ, ψ and Z were changed to remove the time coordinate.
- Writing the information on a file part is updated as other DArray version of the code with println()'s instead of writedlm().

MPIAUX.jl:
- updatePaddings(), transfer_Paddings_f() no longer have "t" as an argument and inside every "t" is removed.

#### v.2.0.0
Version name updated. No other changes. Also this version don't have its own commit. We put here just for book keeping. 

#### v.0.3.8
MPI Code updated --> Z_gl and dZdt_gl are no longer defined. Should improve memory usage.

#### v0.3.6
Parallelization done using MPI with a major change in the structure change of the whole code.
Results match exactly with serial and the julia-parallel version. 

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
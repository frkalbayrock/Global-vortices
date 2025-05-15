# CHANGELOG

#### v0.4.1
We finally added the command line arguments to create the directories where the data is going to be written. Surprisingly, this was always done manually when needed until now because of being extremely lazy.
No other change.

#### v0.4
We tested some new methods on the serial code; mostly about macro usage for fluxes.
We don't see any allocation improvements as long as we @inline the flux functions. (@inline'ing is very important!)
Also have few minor changes to vectorize some calculations and more efficiently. 

All the changes listed:

Evolution.jl:
- Again repositioned some code. 
- Macros for fluxes added in the beginning. (Although, no significant improvement, there is some improvement which we may test in the parallel version later for a bigger lattice to its true effect)
- half_step!() is vecrorized. We also noticed there was an operation (+)  which when it's replaced with (.+) improved the speed and allocations!
- We also have a new leap_forward!() to accomodate the macros used.
- We keep both versions at this stage.

Take-aways to be applied or tested on parallel version:
- Properly vectorize half_step!()
- @inline flux functions if we are gonna keep this version
- Test the macros on parallel version. If we get significant improvements then we use macros, if not we need to definitely use the @inline!


#### v0.3.1
Prepared the code to work with 4-index notation. There was a little problem earlier.(macro implementation didn't work, long story...) This is done to prepare the code for macro implementation which was done earlier but didn't work properly.
Mostly repositioned the code and made sure without changing too much we send in the 4-index "Z" and "dZdt".


#### v0.3
Runge-Kutta applied for learning purposes.

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
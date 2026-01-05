# Global Strings

## Parallel Julia code for 2D+1 lattice simulations for topological defect production with 

#### --->>UNDER CONSTRUCTION<<-----

This numerical simulation code is able to simulate a scattering experiment of two classical wavepackets to create topological defects using the "quantum bridge" that is built into the theory. 


The code uses position-Verlet Leap Frog method for time evolution as the main temporal discretization. However, there are branches that work with Runge-Kutta 4th Order and Iterated Crank-Nicolson with two iterations. 


The novel parallelization method used in the main numerical simulation code is based on high-level Julia language. The different branches leverages MPI and various other distributed computation methods/packages as well. There are also branches with completely new approaches for parallelization that leverages the combination several methods together. For our particular use, we found the fastest method is the one used in the main branch.



This project is a continuation of the following work: 
https://arxiv.org/abs/2308.01962

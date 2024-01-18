This code does vacancy hopping KMC. It supports 3 algorithms:
	- rejection-free KMC
	- first-reaction KMC
	- rejection KMC (currently uses a placeholder accept/reject criterion)

Currently, running a trajectory produces a 3x3 grid plot of the state of the lattice occupancy in 9  X-Y plane cuts.

The user can set various parameters near the beginning of the main section, including:
	Lattice dimensions (Nx, Ny, Nz)
	Temperature (in Kelvin)
	A measure of the bond energy between occupied elements (nn_delta_e)
	KMC algorithm
	Stopping criterion (time or number of steps)

There is also an option to set an "output_time" for how often the 3x3 plot is updated during the KMC trajectory. That variable is set to zero in the current code so every step displays a new image.
 

	

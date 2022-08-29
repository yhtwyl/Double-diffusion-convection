# Double-diffusion-convection
To run this code on linux you have to have gfortran and mpich installed on your machine. These two are easily available online. Once are sure you have these on your system second step would be to run the makefile. In thirs step you have to use >> mpirun -np main_3d.o and link all the remaining subroutines object files. full understanding of the run can be gained from this link https://www.open-mpi.org/doc/current/man1/mpirun.1.php
# I have not tested these codes on Windows or Mac.
# Brief introduction on what these codes simulate
This is a parallelised fortran code for simulating salt fingers in double diffusion convection[1,2]. 
Double diffusive convection (DDC) is the buoyancy-driven flow, with density depending on two different diffusing scalar components, 
distributed such that faster diffusing component gravitationally stabilises the fluid and slower one destabilise it.
The presence of DDC is universal from oceans to stars.
In oceans the two components - heat and salt are convected from top to depths of oceans affecting the climate of earth and life in deep sea.
The stable stratification of molten constituents in the interiors of planets and gases in stars is disturbed due to this instability,
reducing the background gradients to change their internal structure.
This code solves transient Navier-Stokes equations, by following Finite Volume Method and adopting SIMPLER algorithm.
The pnew_3d.inp files contains all the parameters arising from non-dimensionalisatoin of the NS Equation, which can be found in reports[3,4]
In the pnew_3d.inp non-dimensionalised numbers - Rayleigh number, Grashof number, prandtl and schmidt number can be changed to observe the effect of these on the evolution of DDC system. Formulation of the governing non-dimensionalised form of the equation can be found in the MastersThesis.pdf 
The result files are written in binary form and you may need wtite your code to read thees output files.
To expedite the increased size of 3D problem, parallelisation was done using MPI.
At different values of Rayleigh number and density ratio variety of planforms can be observed.  
From the results corresponding to different parameters (supplied in input file 'pnew3d_inp') salt fingers morphology and 
mixing characteristics of double diffusive convection (DDC) for variety of fluids can be observed.

REFERENCES:
1. Huppert, Herbert E., and J. Stewart Turner. "Double-diffusive convection." Journal of Fluid Mechanics 106 (1981): 299-329.
2. Radko, Timour. Double-diffusive convection. Cambridge University Press, 2013.
3. Dhiman, Manoj, Faria Rehman, and O. P. Singh. "Effect of Rayleigh Number and Density Stability Ratio on characteristics of Double-Diffusive Salt Fingers." Journal of Physics: Conference Series. Vol. 759. No. 1. IOP Publishing, 2016.
4. Dhiman, Manoj. Salt fingers in two and three dimensions (MS). Diss. IITMandi, 2016.

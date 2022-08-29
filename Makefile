# makefile: lets try
objects = main_3d.o mpi_3d.o initialise_3d.o exchng3.o fnd3dnbrs.o mpe_decomp1d.o pressure.o vel_coeff3d.o vel_cap3d.o umome_3d.o vmome_3d.o wmome_3d.o mass_3d.o pressure_cor.o uvw_vel.o species_coeff.o species_calc.o tridag.o screen.o duvwdt.o file_name.o

Comp = /usr/local/bin/bin/mpif90 ## You may have to change it depending upon the location of the complier

FFLAGS = -O3 -Wall -Wextra -Wimplicit-interface -fPIC -fmax-errors=1 -g -fcheck=all -fbacktrace ## few of the flags I was interseted in you may change or add new 

main_3d : $(objects)
	$(Comp) -o main_3d $(objects)

%.o : %.f90
	$(Comp) -c $<

#!/bin/bash 
#PBS -q day
#PBS -o out.o
#PBS -e out.e
#PBS -N  pass_time
#PBS -l nodes=1:ppn=8
#PBS -V

cd ${PBS_O_WORKDIR}
echo "Running on: " 
cat ${PBS_NODEFILE} 
cat $PBS_NODEFILE > machines.list
echo "Program Output begins: " 
##/opt/openmpi_intel/bin/mpirun -np 16 -machinefile machines.list 3_main.o
##/opt/soft/share/mpich3p2gcc/bin/mpiexec -np 16 -machinefile machines.list try.o
#/opt/soft/share/openmpigcc/bin/mpiexec -np 16 -machinefile machines.list 3_main.o
mpiexec -np 8 main_3d.o 

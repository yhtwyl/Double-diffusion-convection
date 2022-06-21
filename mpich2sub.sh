#!/bin/bash 
#PBS -j oe 
#PBS -q batch
#PBS -N  pass_time

# it is desirable to request for all processors of
# a node if you have multinode jobs
#PBS -l nodes=1:ppn=8

#PBS -V

cd ${PBS_O_WORKDIR} 
echo "Running on: " 
cat ${PBS_NODEFILE} 
cat ${PBS_NODEFILE}|uniq > node.txt
/opt/mpich2_intel/bin/mpdboot -n 2 -f node.txt 
echo 
echo "Program Output begins: " 
/opt/mpich2_intel/bin/mpiexec -np 8 3_main.o  > output.txt
/opt/mpich2_intel/bin/mpdallexit 

#!/bin/bash 
#PBS -N heat2
##PBS -m abe 
##PBS -M you@northwestern.edu
#PBS -l walltime=00:05:00 
#PBS -q batch 
cd $PBS_O_WORKDIR 
export LD_LIBRARY_PATH=/usr/local/magma/lib:$LD_LIBRARY_PATH
./heat2 heat2.in heat2.out

#!/bin/bash
#PBS -l "walltime=01:30:00"
#PBS -N "fwi.mpi"
#PBS -l "nodes=4:ppn=16:cluster"
#PBS -e job.err
#PBS -o job.out
#PBS -d .
ulimit -c unlimited

date
source environment_sdv.sh
export OMP_NUM_THREADS=16

mpirun -n 2 ./fwi.intel64 ../SetupParams/fwi_params.txt ../SetupParams/fwi_frequencies.txt

date

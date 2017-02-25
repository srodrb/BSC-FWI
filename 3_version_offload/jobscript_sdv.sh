#!/bin/bash
#PBS -l "walltime=00:05:00"
#PBS -N "fwi.mpi"
#PBS -l "nodes=4:ppn=16:cluster"
#PBS -e job.err
#PBS -o job.out
#PBS -d .
ulimit -c unlimited

date
source environment_sdv.sh
export OMP_NUM_THREADS=16

mpiexec.hydra -n 1 ./fwi.intel64 fwi_schedule.txt

date

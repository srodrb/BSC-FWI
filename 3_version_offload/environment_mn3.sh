#!/bin/bash
module purge
module load CMAKE/2.8.12.2
module load gcc/5.1.0
module load intel/15.0.2
module load impi/5.0.1.035
module load ompss/offload

export FWIDIR=$PWD/..
export CC=mpimcc

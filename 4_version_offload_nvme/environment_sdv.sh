#!/bin/bash
module purge
module load cmake/3.6.2
module load gcc/5.3.0-64bit
module load intel/17.0
module load impi/5.1
module load ompss/offload/16.08

export FWIDIR=$PWD/..
export CC=mpimcc

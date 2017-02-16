#!/bin/bash
module purge
module load intel/2017.1
# module load tools/ompss/git
module load cmake 
#INTEL
export CC=icc
#INTEL + OmpSs
# export CC=imcc
export FWIDIR=$PWD/..

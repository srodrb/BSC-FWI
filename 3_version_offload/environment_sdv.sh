#!/bin/bash
module purge
module load CMAKE
module load gcc
module load intel
module load impi
module load ompss/offload

export FWIDIR=$PWD/..
export CC=mpimcc

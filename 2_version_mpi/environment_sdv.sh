#!/bin/bash
module purge
module load intel
module load impi
module load CMAKE

export CC=mpiicc
export FWIDIR=$PWD/..

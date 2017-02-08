#!/bin/bash
module purge
module load intel/16.0.0
module load impi/5.1.3.210
module load CMAKE/2.8.12.2

export CC=mpiicc
export FWIDIR=$PWD/..

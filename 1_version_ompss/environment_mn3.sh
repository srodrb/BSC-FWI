#!/bin/bash
module purge
module load gcc/6.1.0
module load intel/2017.1
module load tools/ompss/git
module load CMAKE/2.8.12.2

export CC=imcc
export FWIDIR=$PWD/..

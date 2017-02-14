#!/bin/bash
module purge
module load intel/2017.1
module load cmake 

export CC=icc
export FWIDIR=$PWD/..

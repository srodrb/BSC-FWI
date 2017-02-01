#!/bin/bash
module purge
module load intel/16.0.0
module load CMAKE/2.8.12.2

export CC=icc
export FWIDIR=$PWD

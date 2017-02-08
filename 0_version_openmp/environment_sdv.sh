#!/bin/bash
module purge
module load intel
module load CMAKE

export CC=icc
export FWIDIR=$PWD/..

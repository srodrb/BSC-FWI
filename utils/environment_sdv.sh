#!/bin/bash
module purge
module load cmake/3.6.2
module load gcc
export CC=gcc
export FWIDIR=$PWD/..

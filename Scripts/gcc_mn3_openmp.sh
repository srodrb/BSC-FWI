#!/bin/bash
module purge
module load gcc/6.1.0
module load CMAKE/2.8.12.2

export CC=gcc
export FWIDIR=$PWD/..

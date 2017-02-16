#!/bin/bash
#SBATCH -J mockup-fwi
#SBATCH --qos=gen
#SBATCH --time=01:00:00
#SBATCH -D .
#SBATCH -e job.knl.err
#SBATCH -o job.knl.out
#SBATCH --ntasks=1 
#SBATCH -N 1
#SBATCH --exclusive
#SBATCH --nodelist=knl04
#SBATCH --constraint=cache

sleep 1h


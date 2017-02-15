#!/usr/bin
export FWIDIR=$PWD/..
export OMP_NUM_THREADS=16
export KMP_AFFINITY=verbose,compact
# export KMP_AFFINITY=granularity=fine,compact,verbose
rm fwi.log
# gdb --args ./fwi.bin fwi_schedule.txt
numactl --membind=0,1 taskset -c 0-15 ./fwi.bin fwi_schedule.test
# ./fwi.bin fwi_schedule.test

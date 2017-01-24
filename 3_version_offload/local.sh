#!/bin/bash
# Resource requirements
# SN : number of slave nodes
# PPN_SN : number of processes per slave node
# WN : number of worker nodes
# PPN_WN: number of processes per worker node
# Total amount of resources: max(sn) + max(wn) + 1


export SN=1
export PPN_SN=1
	export WN=3
export PPN_WN=1

mkdir -p .tmpdir

export NX_OFFL_HOSTFILE=".tmpdir/slave.hosts"
export OFFL_NX_OFFL_HOSTFILE=".tmpdir/worker.hosts"
rm -f $NX_OFFL_HOSTFILE
for h in $(seq $SN); do
	  $( echo localhost >> $NX_OFFL_HOSTFILE )
done

rm -f $OFFL_NX_OFFL_HOSTFILE
for h in $(seq $WN); do
  $( echo localhost >> $OFFL_NX_OFFL_HOSTFILE )
done

export NX_SMP_WORKERS=1
export OFFL_NX_SMP_WORKERS=1

# Unset BOOTSTRAP=lsf. LSF does not support dynamic process execution completely.
bin=fwi.intel64
#export TERMINAL="xterm -e"
#export TERMINAL=
#export GDB="gdb --args"
#export GDB="valgrind --tool=memcheck"
#export NX_OFFL_GDB=1
#export NX_OFFL_DEBUG=5

#export MPICH_DBG_CLASS="PT2PT"
#export MPICH_DBG_LEVEL="VERBOSE"

#MPIRUN=mpiexec.hydra
#MPIRUN="mpiexec.hydra -outfile-pattern=output-%g-%r -errfile-pattern=output-%g-%p-%r"
MPIRUN="mpiexec.hydra -prepend-pattern=[%g-%r]"

export OMP_STACKSIZE="3M"
#export NX_ARGS="--disable-ut"

${MPIRUN} -genv -wdir $(pwd) -np 1 $TERMINAL $GDB ./$bin fwi_schedule.debug


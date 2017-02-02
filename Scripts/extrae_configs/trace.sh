#!/bin/bash

source /apps/BSCTOOLS/extrae/3.4.1/impi_mt+libgomp4.2/etc/extrae.sh

export EXTRAE_CONFIG_FILE=extrae.xml
export NX_ARGS="${NX_ARGS} --instrumentation=extrae "
export LD_PRELOAD=${EXTRAE_HOME}/lib/libnanosmpitrace.so
#export LD_PRELOAD=${EXTRAE_HOME}/lib/libnanosmpitracef.so

$*

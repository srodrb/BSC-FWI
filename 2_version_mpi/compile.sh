#!/bin/bash
mpiicc -std=c99 -restrict -qopenmp  -g -Wall -DDISTRIBUTED_MEMORY_IMPLEMENTATION -I../common -I. ../common/fwi_constants.c ../common/fwi_sched.c ../common/fwi_common.c ../common/fwi_kernel.c ./fwi_propagator.c ./fwi_main.c -o fwi.intel64 -lm

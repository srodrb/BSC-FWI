#!/bin/bash
mpimcc -DDEBUG -DDISTRIBUTED_MEMORY_IMPLEMENTATION -g -Wall --ompss --no-copy-deps -I../common -I. ../common/fwi_constants.c ../common/fwi_sched.c ../common/fwi_common.c ../common/fwi_kernel.c ./fwi_offload.c ./fwi_propagator.c ./fwi_main.c -o fwi.intel64 -lm

# possible values <YES,NO>
PERFORM_IO ?=NO
# possible values <YES,NO>
IO_STATS   ?=NO
# possible values <YES,NO>
DEBUG      ?=NO
# Uset local NVME local IO? <YES,NO>
USE_NVME   ?=NO


ifeq ($(PERFORM_IO),NO)
	DEFINES  += -DDO_NOT_PERFORM_IO
endif

ifeq ($(IO_STATS),NO)
	DEFINES  +=
	CFLAGS   += -Wno-unused-variable -Wno-unused-but-set-variable
else
	DEFINES  += -DLOG_IO_STATS
	CFLAGS   +=
endif

ifeq ($(DEBUG),YES)
	DEFINES  += -DDEBUG
	CFLAGS   += -g
endif

ifeq ($(USE_NMVE),YES)
	DEFINES  += -DUSE_NMVE
endif




# Compiler and libraries setup
MPICC    =mpimcc --Wn,-xHost
CC       =icc
CFLAGS   =-Iinclude -DMPI_ENABLED -O3 -std=gnu99
MCFLAGS  =--ompss --no-copy-deps
LDFLAGS  =-lrt -lm

TARGET=fwi.intel64 create_schedule ModelGenerator

all:$(TARGET)

fwi.intel64:fwi_main.o fwi_common.o fwi_kernel.o fwi_propagator.o fwi_offload.o fwi_sched.o
	$(MPICC) $(DEFINES) $(CFLAGS) $(MCFLAGS) $(LDFLAGS) $^ -o $@

fwi_main.o:fwi_main.c
	$(MPICC) $(DEFINES) $(CFLAGS) $(MCFLAGS) -c $^ -o $@

fwi_common.o:fwi_common.c
	$(MPICC) $(DEFINES) $(CFLAGS) $(MCFLAGS) -c $^ -o $@

fwi_kernel.o:fwi_kernel.c
	$(MPICC) $(DEFINES) $(CFLAGS) $(MCFLAGS) -c $^ -o $@

fwi_propagator.o:fwi_propagator.c
	$(MPICC) $(DEFINES) $(CFLAGS) $(MCFLAGS) -c $^ -o $@

fwi_offload.o:fwi_offload.c
	$(MPICC) $(DEFINES) $(CFLAGS) $(MCFLAGS) -c $^ -o $@

fwi_sched.o:fwi_sched.c
	$(MPICC) $(DEFINES) $(CFLAGS) $(MCFLAGS) -c $^ -o $@

create_schedule:fwi_schedule.c
	$(CC) $(DEFINES) -std=gnu99 fwi_schedule.c -o create_schedule -lm

ModelGenerator:fwi_generatemodel.c fwi_common.o
	$(MPICC) $(DEFINES) $(CFLAGS) $(MCFLAGS) $(LDFLAGS) fwi_generatemodel.c fwi_common.o -o ModelGenerator -lm

.PHONY:all clean input submit schedule cleanout

clean:
	rm -rf *.o *.optrpt $(TARGET) fname* InputModels/* results/*

input:
	./ModelGenerator ../SetupParams/fwi_params.txt ../SetupParams/fwi_frequencies.txt

submit:
	qsub < sdv_x2x.pbs && watch qstat

schedule:
	./create_schedule fwi_params.txt fwi_frequencies.txt


cleanout:
	rm -rf worker.* slave.* *.err *.out core*

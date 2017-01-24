# possible values <INTEL,GNU>
COMPILER   ?=INTEL
# possible values <YES,NO>
PERFORM_IO ?=YES
# possible values <YES,NO>
IO_STATS   ?=YES
# possible values <YES,NO>
DEBUG      ?=NO
# Num MPI ranks 
NCPUS      ?= 2

ifeq ($(COMPILER),INTEL)
	CC        = icc
	MPICC     = mpiicc
	MPICFLAGS = -std=c99 -xHost -qopenmp -restrict -prec-div -fp-model precise -ftz -qopt-prefetch-distance=64
	# MPICFLAGS = -std=c99 -O2 -qopenmp -restrict
	CFLAGS    = -O0 -std=c99 -qopenmp -restrict -DNO_MPI
	CINCS     =
	CLIBS     = -lm
else
	CC        = gcc
	MPICC     = mpicc
	MPICFLAGS = -std=c99 -O3 -march=native -fopenmp -Wall #-fopt-info -fopt-info-vec-missed
	CFLAGS    = -O0 -std=c99 -qopenmp -restrict -DNO_MPI
	CINCS     =
	CLIBS     = -lm
endif

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



# Directories, source, include, Objects, etc..
MPIBUILDDIR = mpibuild
SEQBUILDDIR = seqbuild

TARGET  = fwi.intel64


all:$(TARGET) ModelGenerator ModelGenerator ScheduleGenerator

fwi.intel64:mpi_fwi_main.o mpi_fwi_common.o mpi_fwi_kernel.o mpi_fwi_propagator.o
	$(MPICC) $(DEFINES) $(MPICFLAGS) $(CINCS) $+ -o $@ $(CLIBS)

mpi_fwi_main.o: fwi_main.c
	$(MPICC) $(DEFINES) $(MPICFLAGS) $(CINCS) -c $+ -o $@

mpi_fwi_common.o: fwi_common.c
	$(MPICC) $(DEFINES) $(MPICFLAGS) $(CINCS) -c $+ -o $@

mpi_fwi_kernel.o: fwi_kernel.c
	$(MPICC) $(DEFINES) $(MPICFLAGS) $(CINCS) -c $+ -o $@

mpi_fwi_propagator.o: fwi_propagator.c
	$(MPICC) $(DEFINES) $(MPICFLAGS) $(CINCS) -c $+ -o $@


# compilation of tools (do not use MPI)
ModelGenerator:fwi_generatemodel.o seq_fwi_common.o seq_fwi_kernel.o seq_fwi_propagator.o
	$(CC) $(DEFINES) $(CFLAGS) $(CINCS) $+ -o $@ $(CLIBS)

ScheduleGenerator:fwi_schedule.o
	$(CC) $(DEFINES) $(CFLAGS) $(CINCS) $+ -o $@ $(CLIBS)

fwi_schedule.o:fwi_schedule.c
	$(CC) $(DEFINES) $(CFLAGS) $(CINCS) -c $+ -o $@

fwi_generatemodel.o:fwi_generatemodel.c
	$(CC) $(DEFINES) $(CFLAGS) $(CINCS) -c $+ -o $@

seq_fwi_common.o: fwi_common.c
	$(CC) $(DEFINES) $(CFLAGS) $(CINCS) -c $+ -o $@

seq_fwi_kernel.o: fwi_kernel.c
	$(CC) $(DEFINES) $(CFLAGS) $(CINCS) -c $+ -o $@

seq_fwi_propagator.o: fwi_propagator.c
	$(CC) $(DEFINES) $(CFLAGS) $(CINCS) -c $+ -o $@


.PHONY:all clean run irun debug memcheck input

clean:
	rm -rf *.o *.optrpt *.log $(TARGET) *.row *.prv *.pcf InputModels/* Results/* *.err *.out core* mpi_*

input: ModelGenerator
	./ModelGenerator ../SetupParams/fwi_params.txt ../SetupParams/fwi_frequencies.txt

schedule:
	./ScheduleGenerator fwi_params.txt fwi_frequencies.txt

run: fwi.intel64
	bsub < jobscript.sh

irun: fwi.intel64
	mpirun -np $(NCPUS) ./fwi.intel64 ../SetupParams/fwi_params.txt ../SetupParams/fwi_frequencies.txt


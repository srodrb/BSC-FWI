# possible values <INTEL,GNU>
BACKEND   ?=INTEL
# possible values <YES,NO>
PERFORM_IO ?=NO
# possible values <YES,NO>
IO_STATS   ?=NO
# possible values <YES,NO>
DEBUG      ?=NO

# Add include directory
CINCS =-I../common -I./

ifeq ($(BACKEND),INTEL)
	CC        = imcc --ompss -O3
	# CFLAGS    = -std=c99 -O2 -restrict -qopenmp #-qopt-report=2 -opt-streaming-cache-evict=0 -opt-streaming-stores always
	# CFLAGS    = --Wn,-march=native --Wn,-restrict --Wn,-prec-div --Wn,-fp-model=precise --Wn,-ftz --Wn,-qopt-prefetch-distance=64
	CFLAGS    = --Wn,-xHost --Wn,-restrict
	CINCS     +=
	CLIBS     = -lm
else
	CC        = mcc --ompss -O3
	CFLAGS    = --Wn,-march=native --Wn,-restrict
	CINCS     +=
	CLIBS     = -lm
endif

ifeq ($(PERFORM_IO),NO)
	DEFINES  += -DDO_NOT_PERFORM_IO
endif

ifeq ($(IO_STATS),NO)
	DEFINES  +=
	CFLAGS   += --Wn,-Wno-unused-variable --Wn,-Wno-unused-but-set-variable
else
	DEFINES  += -DLOG_IO_STATS
	CFLAGS   +=
endif

ifeq ($(DEBUG),YES)
	DEFINES  += -DDEBUG
	CFLAGS   += -g
endif

TARGET=fwi.intel64 ModelGenerator

all:$(TARGET)

fwi.intel64:fwi_main.o fwi_common.o fwi_kernel.o fwi_propagator.o
	$(CC) $(DEFINES) $(CFLAGS) $(CINCS) $+ -o $@ $(CLIBS)

ModelGenerator:fwi_generatemodel.o fwi_common.o fwi_kernel.o fwi_propagator.o
	$(CC) $(DEFINES) $(CFLAGS) $(CINCS) $+ -o $@ $(CLIBS)

fwi_generatemodel.o:../common/fwi_generatemodel.c
	$(CC) $(DEFINES) $(CFLAGS) $(CINCS) $^ -c -o $@

fwi_main.o:fwi_main.c
	$(CC) $(DEFINES) $(CFLAGS) $(CINCS) $^ -c -o $@

fwi_common.o:../common/fwi_common.c
	$(CC) $(DEFINES) $(CFLAGS) $(CINCS) $^ -c -o $@

fwi_kernel.o:../common/fwi_kernel.c
	$(CC) $(DEFINES) $(CFLAGS) $(CINCS) $^ -c -o $@

fwi_propagator.o:fwi_propagator.c
	$(CC) $(DEFINES) $(CFLAGS) $(CINCS) $^ -c -o $@


.PHONY:all clean run debug memcheck input printvars

clean:
	rm -rf *.o *.optrpt *.log $(TARGET) *.log *.csv ../Results/*

input: ModelGenerator
	./ModelGenerator ../SetupParams/fwi_params.txt ../SetupParams/fwi_frequencies.txt

run: fwi.intel64
	./fwi.intel64 ../SetupParams/fwi_params.txt ../SetupParams/fwi_frequencies.txt

debug: fwi.intel64
	gdb --args fwi.intel64 ../SetupParams/fwi_params.txt ../SetupParams/fwi_frequencies.txt

memcheck: fwi.intel64
	valgrind --leak-check=full ./fwi.intel64 ../SetupParams/fwi_params.txt ../SetupParams/fwi_frequencies.txt > valgrind.log 2>&1


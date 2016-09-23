# possible values <INTEL,GNU>
COMPILER   ?=GNU
# possible values <YES,NO>
PERFORM_IO ?=YES
# possible values <YES,NO>
IO_STATS   ?=NO
# possible values <YES,NO>
DEBUG      ?=NO
# Num MPI ranks 
NCPUS ?= 2

ifeq ($(COMPILER),INTEL)
	CC        = icc
	MPICC     = mpicc
	# CFLAGS    = -std=c99 -O2 -restrict -qopenmp #-qopt-report=2
	# -opt-streaming-cache-evict=0 -opt-streaming-stores always
	CFLAGS    = -std=c99 -xHost -qopenmp -restrict -prec-div -fp-model precise -ftz -opt-prefetch-distance=64
	CINCS     =
	CLIBS     = -lm
else
	CC        = gcc
	MPICC     = mpicc
	CFLAGS    = -std=c99 -O3 -march=native -fopenmp -Wall #-fopt-info -fopt-info-vec-missed
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

TARGET=fwi.intel64 ModelGenerator

all:$(TARGET)

fwi.intel64:fwi_main.o fwi_common.o fwi_kernel.o fwi_propagator.o
	$(MPICC) $(DEFINES) $(CFLAGS) $(CINCS) $+ -o $@ $(CLIBS) 

ModelGenerator:fwi_generatemodel.o fwi_common.o fwi_kernel.o fwi_propagator.o
	$(MPICC) $(DEFINES) $(CFLAGS) $(CINCS) $+ -o $@ $(CLIBS)

fwi_generatemodel.o:fwi_generatemodel.c 
	$(MPICC) $(DEFINES) $(CFLAGS) $(CINCS) $+ -c

fwi_main.o:fwi_main.c 
	$(MPICC) $(DEFINES) $(CFLAGS) $(CINCS) $+ -c

fwi_common.o:fwi_common.c
	$(MPICC) $(DEFINES) $(CFLAGS) $(CINCS) $+ -c

fwi_kernel.o:fwi_kernel.c
	$(MPICC) $(DEFINES) $(CFLAGS) $(CINCS) $+ -c

fwi_propagator.o:fwi_propagator.c
	$(MPICC) $(DEFINES) $(CFLAGS) $(CINCS) $+ -c


.PHONY:all clean run irun debug memcheck input

clean:
	rm -rf *.o *.optrpt *.log $(TARGET) *.row *.prv *.pcf InputModels/* Results/* *.err *.out core* mpi_*

input: ModelGenerator
	./ModelGenerator ../SetupParams/fwi_params.txt ../SetupParams/fwi_frequencies.txt

run: fwi.intel64
	bsub < jobscript.sh

irun: fwi.intel64
	mpirun -np $(NCPUS) ./fwi.intel64 ../SetupParams/fwi_params.txt ../SetupParams/fwi_frequencies.txt

debug: fwi.intel64
	mpirun -np $(NCPUS) xterm -e gdb --args fwi.intel64 ../SetupParams/fwi_params.txt ../SetupParams/fwi_frequencies.txt

memcheck: fwi.intel64
	valgrind ./fwi.intel64 ../SetupParams/fwi_params.txt ../SetupParams/fwi_frequencies.txt

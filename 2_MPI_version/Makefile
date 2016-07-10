CC      = gcc
MPICC   = mpicc
CFLAGS  = -std=c99 -O3
INCS    =
LIBS    = -lm
DEFINES = -DDO_NOT_PERFOM_IO -DDEBUG -DDEBUG_MPI -DDEBUG_TO_STDERR
LDFLAGS = 

ifeq ($(CC), gcc)
    CFLAGS += -fopenmp -march=native -fopt-info #-fopt-info-vec-missed
else
ifeq ($(CC), icc)
    CFLAGS += -openmp -mavx -qopt-report=2
   #CFLAGS += -openmp -xHost -qopt-assume-safe-padding
else
ifeq ($(CC), mcc)
    CFLAGS = --opmss --Wn,-std=gnu99 --Wn,-xHost
endif
endif
endif

NCPUS ?= 2

TARGET=fwi.intel64 ModelGenerator

all:$(TARGET)

fwi.intel64:fwi_main.o fwi_common.o fwi_kernel.o fwi_propagator.o
	$(MPICC) $(DEFINES) $(CFLAGS) $(INCS) $+ -o $@ $(LIBS) $(LDFLAGS)

ModelGenerator:fwi_generatemodel.o fwi_common.o fwi_kernel.o fwi_propagator.o
	$(MPICC) $(DEFINES) $(CFLAGS) $(INCS) $+ -o $@ $(LIBS) $(LDFLAGS)

fwi_generatemodel.o:fwi_generatemodel.c 
	$(MPICC) $(DEFINES) $(CFLAGS) $(INCS) $+ -c

fwi_main.o:fwi_main.c 
	$(MPICC) $(DEFINES) $(CFLAGS) $(INCS) $+ -c

fwi_common.o:fwi_common.c
	$(MPICC) $(DEFINES) $(CFLAGS) $(INCS) $+ -c

fwi_kernel.o:fwi_kernel.c
	$(MPICC) $(DEFINES) $(CFLAGS) $(INCS) $+ -c

fwi_propagator.o:fwi_propagator.c
	$(MPICC) $(DEFINES) $(CFLAGS) $(INCS) $+ -c


.PHONY:all clean run debug memcheck input

clean:
	rm -rf *.o *.optrpt $(TARGET) *.row *.prv *.pcf InputModels/* results/* *.err *.out core* mpi_*

input: ModelGenerator
	./ModelGenerator ../SetupParams/fwi_params.txt ../SetupParams/fwi_frequencies.txt

run: fwi.intel64
	bsub < jobscript.sh

irun: fwi.intel64
	mpirun -np $(NCPUS) ./fwi.intel64 ../SetupParams/fwi_params.txt ../SetupParams/fwi_frequencies.txt


debug: fwi.intel64
	gdb --args fwi.intel64 ../SetupParams/fwi_params.txt ../SetupParams/fwi_frequencies.txt

memcheck: fwi.intel64
	valgrind ./fwi.intel64 ../SetupParams/fwi_params.txt ../SetupParams/fwi_frequencies.txt

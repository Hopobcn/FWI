# possible values <PGI>
COMPILER    =PGI
# possible values <YES,NO>
PERFORM_IO ?=YES
# possible values <YES,NO>
IO_STATS   ?=NO
# possible values <YES,NO>
DEBUG      ?=NO
# possible values <YES,NO>
USE_CUDA   ?=YES
# possible values <YES,NO>
TRACE_CUDA ?=NO
# Num MPI ranks 
NCPUS ?= 2

HOST=$(shell hostname)

ifeq ($(COMPILER),PGI)
	CC        = pgcc
	NVCC      = nvcc
	MPICC     = mpicc
	
	ifeq ($(HOST),hulk)
		HOST_CC   = gcc-4.9
	else
		HOST_CC   = /apps/GCC/4.9.3/bin/gcc
	endif
	
	CFLAGS    =  -c99 -fast -mp -acc -ta=tesla,cc20,cc35,cc50,maxregcount:128,lineinfo -Minfo=accel,inline,mp -Minline=IDX
	NVFLAGS   = -std=c++11 -ccbin $(HOST_CC) -O3 -gencode arch=compute_30,code=\"sm_35,sm_37\" -gencode arch=compute_50,code=sm_52 -res-usage -lineinfo
	CINCS     = -I$(CUDA_HOME)/include -I$(MPI_HOME)/include $(shell pkg-config --cflags hwloc) #-I$(ACC_HOME)/include
	CLIBS     = -lm -L$(CUDA_HOME)/lib64 -lcuda -lcudart -L$(MPI_HOME)/lib64 -lmpi $(shell pkg-config --libs hwloc) -lstdc++
else
	$(error $(COMPILER) not supported)
endif

ifeq ($(PERFORM_IO),NO)
	DEFINES  += -DDO_NOT_PERFORM_IO
endif

ifeq ($(IO_STATS),NO)
	DEFINES  +=
	CFLAGS   += #-Wno-unused-variable -Wno-unused-but-set-variable
else
	DEFINES  += -DLOG_IO_STATS
	CFLAGS   +=
endif

ifeq ($(DEBUG),YES)
	DEFINES  += -DDEBUG
	CFLAGS   += -g
endif

ifeq ($(USE_CUDA),YES)
	DEFINES += -DUSE_CUDA

	# SCELL optimized kernels:
	DEFINES += -DSCELL_TL
	DEFINES += -DSCELL_BL
	#DEFINES += -DSCELL_BR
	#DEFINES += -DSCELL_TR

	# VCELL optimized kernels:
	DEFINES += -DVCELL_TL
endif

ifeq ($(TRACE_CUDA),YES)
	DEFINES += -DTRACE_CUDA
	CLIBS   += -lnvToolsExt
endif

TARGET=fwi.intel64 ModelGenerator

all:$(TARGET)

fwi.intel64: fwi_main.o fwi_common.o fwi_kernel.o fwi_propagator.o fwi_common_cuda.o fwi_propagator_cuda.o
	$(MPICC) $(DEFINES) $(CFLAGS) $(CINCS) $+ -o $@ $(CLIBS)

ModelGenerator: fwi_generatemodel.o fwi_common.o fwi_kernel.o fwi_propagator.o fwi_common_cuda.o fwi_propagator_cuda.o
	$(MPICC) $(DEFINES) $(CFLAGS) $(CINCS) $+ -o $@ $(CLIBS)

fwi_generatemodel.o: fwi_generatemodel.c
	$(MPICC) $(DEFINES) $(CFLAGS) $(CINCS) $+ -c

fwi_main.o: fwi_main.c
	$(MPICC) $(DEFINES) $(CFLAGS) $(CINCS) $+ -c

fwi_common.o: fwi_common.c
	$(MPICC) $(DEFINES) $(CFLAGS) $(CINCS) $+ -c

fwi_kernel.o: fwi_kernel.c
	$(MPICC) $(DEFINES) $(CFLAGS) $(CINCS) $+ -c

fwi_propagator.o: fwi_propagator.c
	$(MPICC) $(DEFINES) $(CFLAGS) $(CINCS) $+ -c

fwi_common_cuda.o: fwi_common.cu
	$(NVCC) $(DEFINES) $(NVFLAGS) $(CINCS) $^ -c -o $@ 2>&1| c++filt --no-params

fwi_propagator_cuda.o: fwi_propagator.cu
	$(NVCC) $(DEFINES) $(NVFLAGS) $(CINCS) $^ -c -o $@ 2>&1| c++filt --no-params


.PHONY:all clean run irun debug memcheck input

clean:
	rm -rf *.o *.optrpt *.log $(TARGET) *.row *.prv *.pcf InputModels/* results/* *.err *.out core* mpi_*

input: ModelGenerator
	./ModelGenerator ../SetupParams/fwi_params.txt ../SetupParams/fwi_frequencies.txt

run: fwi.intel64
	bsub < jobscript_mt.sh

irun: fwi.intel64
	mpirun -np $(NCPUS) ./fwi.intel64 ../SetupParams/fwi_params.txt ../SetupParams/fwi_frequencies.txt

debug: fwi.intel64
	mpirun -np $(NCPUS) xterm -e cuda-gdb --args fwi.intel64 ../SetupParams/fwi_params.txt ../SetupParams/fwi_frequencies.txt

memcheck: fwi.intel64
	valgrind ./fwi.intel64 ../SetupParams/fwi_params.txt ../SetupParams/fwi_frequencies.txt
profile: profile-gpu

profile-cpu: fwi.intel64
	mpirun -np $(NCPUS) nvprof -o fwi.%q{OMPI_COMM_WORLD_RANK}.nvprof --cpu-profiling on --cpu-profiling-mode top-down ./fwi.intel64 ../SetupParams/fwi_params.txt ../SetupParams/fwi_frequencies.txt

profile-gpu: fwi.intel64
	mpirun -np $(NCPUS) nvprof -o fwi.%q{OMPI_COMM_WORLD_RANK}.nvprof ./fwi.intel64 ../SetupParams/fwi_params.txt ../SetupParams/fwi_frequencies.txt


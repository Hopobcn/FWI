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

ifeq ($(COMPILER),PGI)
	CC        = pgcc
	NVCC      = nvcc
	CFLAGS    = -c99 -fast -mp -acc -ta=tesla,cc20,cc35,cc50,maxregcount:128,lineinfo -Minfo=accel,inline,mp -Minline=IDX
	NVFLAGS   = -std=c++11 -O3 -gencode arch=compute_20,code=sm_20 -gencode arch=compute_30,code=\"sm_35,sm_37\" -gencode arch=compute_50,code=sm_52 -res-usage -lineinfo
	CINCS     = -I$(CUDA_HOME)/include #-I$(ACC_HOME)/include
	CLIBS     = -lm -L$(CUDA_HOME)/lib64 -lcuda -lcudart
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
endif

ifeq ($(TRACE_CUDA),YES)
	DEFINES += -DTRACE_CUDA
	CLIBS   += -lnvToolsExt
endif

TARGET=fwi.intel64 ModelGenerator

all:$(TARGET)

fwi.intel64: fwi_main.o fwi_common.o fwi_kernel.o fwi_propagator.o fwi_propagator_cuda.o
	$(CC) $(DEFINES) $(CFLAGS) $(CINCS) $+ -o $@ $(CLIBS)

ModelGenerator: fwi_generatemodel.o fwi_common.o fwi_kernel.o fwi_propagator.o fwi_propagator_cuda.o
	$(CC) $(DEFINES) $(CFLAGS) $(CINCS) $+ -o $@ $(CLIBS)

fwi_generatemodel.o: fwi_generatemodel.c
	$(CC) $(DEFINES) $(CFLAGS) $(CINCS) $^ -c -o $@

fwi_main.o: fwi_main.c
	$(CC) $(DEFINES) $(CFLAGS) $(CINCS) $^ -c -o $@

fwi_common.o: fwi_common.c
	$(CC) $(DEFINES) $(CFLAGS) $(CINCS) $^ -c -o $@

fwi_kernel.o: fwi_kernel.c
	$(CC) $(DEFINES) $(CFLAGS) $(CINCS) $^ -c -o $@

fwi_propagator.o: fwi_propagator.c
	$(CC) $(DEFINES) $(CFLAGS) $(CINCS) $^ -c -o $@

fwi_propagator_cuda.o: fwi_propagator.cu
	$(NVCC) $(DEFINES) $(NVFLAGS) $(CINCS) $^ -c -o $@ 2>&1| c++filt --no-params


.PHONY:all clean run debug memcheck input

clean:
	rm -rf *.o *.optrpt *.log $(TARGET) fname* InputModels/* results/*

input: ModelGenerator
	./ModelGenerator ../SetupParams/fwi_params.txt ../SetupParams/fwi_frequencies.txt

run: fwi.intel64
	./fwi.intel64 ../SetupParams/fwi_params.txt ../SetupParams/fwi_frequencies.txt

debug: fwi.intel64
	cuda-gdb --args fwi.intel64 ../SetupParams/fwi_params.txt ../SetupParams/fwi_frequencies.txt

memcheck: fwi.intel64
	valgrind --leak-check=full ./fwi.intel64 ../SetupParams/fwi_params.txt ../SetupParams/fwi_frequencies.txt > valgrind.log 2>&1

profile: profile-gpu

profile-cpu: fwi.intel64
	nvprof --cpu-profiling on --cpu-profiling-mode top-down ./fwi.intel64 ../SetupParams/fwi_params.txt ../SetupParams/fwi_frequencies.txt

profile-gpu: fwi.intel64
	nvprof ./fwi.intel64 ../SetupParams/fwi_params.txt ../SetupParams/fwi_frequencies.txt


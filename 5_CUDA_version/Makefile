CC       = pgcc
NVCC     = nvcc
CFLAGS   = -g -c99 -fast -mp -acc -ta=tesla,cc20,cc35,cc50,maxregcount:128,lineinfo -Minfo=accel,inline,mp -Minline=IDX -Mcuda
NVFLAGS  = -g -std=c++11 -O3 -gencode arch=compute_20,code=sm_20 -gencode arch=compute_30,code=\"sm_35,sm_37\" -gencode arch=compute_50,code=sm_52 -res-usage -lineinfo
INCS     = -I$(CUDA_HOME)/include #-I$(ACC_HOME)/include
LIBS     = -lm
DEFINES  += -DDO_NOT_PERFORM_IO -DUSE_CUDA







TARGET=fwi.intel64 ModelGenerator

all:$(TARGET)

fwi.intel64: fwi_main.o fwi_common.o fwi_kernel.o fwi_propagator.o fwi_propagator_cuda.o
	$(CC) $(DEFINES) $(CFLAGS) $(INCS) $+ -o $@ $(LIBS)

ModelGenerator: fwi_generatemodel.o fwi_common.o fwi_kernel.o fwi_propagator.o fwi_propagator_cuda.o
	$(CC) $(DEFINES) $(CFLAGS) $(INCS) $+ -o $@ $(LIBS)

fwi_generatemodel.o: fwi_generatemodel.c
	$(CC) $(DEFINES) $(CFLAGS) $(INCS) $^ -c -o $@

fwi_main.o: fwi_main.c
	$(CC) $(DEFINES) $(CFLAGS) $(INCS) $^ -c -o $@

fwi_common.o: fwi_common.c
	$(CC) $(DEFINES) $(CFLAGS) $(INCS) $^ -c -o $@

fwi_kernel.o: fwi_kernel.c
	$(CC) $(DEFINES) $(CFLAGS) $(INCS) $^ -c -o $@

fwi_propagator.o: fwi_propagator.c
	$(CC) $(DEFINES) $(CFLAGS) $(INCS) $^ -c -o $@

fwi_propagator_cuda.o: fwi_propagator.cu
	$(NVCC) $(DEFINES) $(NVFLAGS) $(INCS) $^ -c -o $@ 2>&1| c++filt --no-params


.PHONY:all clean run debug memcheck input

clean:
	rm -rf *.o *.optrpt $(TARGET) fname* InputModels/* results/*

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


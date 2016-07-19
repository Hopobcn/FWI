CC       = gcc
CFLAGS   = -std=c99 -g -O3
INCS     =
LIBS     = -lm
DEFINES  += #-DDEBUG

ifeq ($(CC), gcc)
    CFLAGS += -fopenmp -march=native -Wall # -fopt-info #-Wno-unused-result
else
ifeq ($(CC), icc)
    CFLAGS += -openmp -mavx -qopt-report=2
   #CFLAGS += -openmp -xHost -qopt-report=2
endif
endif

TARGET=fwi.test

all:$(TARGET)

fwi.test:fwi_test.o fwi_common.o fwi_kernel.o #fwi_propagator.o
	$(CC) $(DEFINES) $(CFLAGS) $(INCS) $+ -o $@ $(LIBS)

#ModelGenerator:fwi_generatemodel.o fwi_common.o fwi_kernel.o fwi_propagator.o
#	$(CC) $(DEFINES) $(CFLAGS) $(INCS) $+ -o $@ $(LIBS)

#fwi_generatemodel.o:fwi_generatemodel.c
#	$(CC) $(DEFINES) $(CFLAGS) $(INCS) $^ -c -o $@

fwi_test.o:fwi_test.c
	$(CC) $(DEFINES) $(CFLAGS) $(INCS) $^ -c -o $@

fwi_common.o:fwi_common.c
	$(CC) $(DEFINES) $(CFLAGS) $(INCS) $^ -c -o $@

fwi_kernel.o:fwi_kernel.c
	$(CC) $(DEFINES) $(CFLAGS) $(INCS) $^ -c -o $@

#fwi_propagator.o:fwi_propagator.c
#	$(CC) $(DEFINES) $(CFLAGS) $(INCS) $^ -c -o $@


.PHONY:all clean test

clean:
	rm -rf *.o *.optrpt $(TARGET) fname* InputModels/* results/*

test: fwi.test
	./fwi.test ../SetupParams/fwi_params.txt ../SetupParams/fwi_frequencies.txt ../0_OpenMP_version ../5_CUDA_version


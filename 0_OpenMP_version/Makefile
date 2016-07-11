CC       = gcc
CFLAGS   = -std=c99 -g -O3
INCS     =
LIBS     = -lm
DEFINES  += -DDO_NOT_PERFORM_IO #-DDEBUG

ifeq ($(CC), gcc)
    CFLAGS += -fopenmp -march=native -Wall -fopt-info #-Wno-unused-result
else
ifeq ($(CC), icc)
    CFLAGS += -openmp -mavx -qopt-report=2
   #CFLAGS += -openmp -xHost -qopt-report=2
endif
endif

TARGET=fwi.intel64 ModelGenerator

all:$(TARGET)

fwi.intel64:fwi_main.o fwi_common.o fwi_kernel.o fwi_propagator.o
	$(CC) $(DEFINES) $(CFLAGS) $(INCS) $+ -o $@ $(LIBS)

ModelGenerator:fwi_generatemodel.o fwi_common.o fwi_kernel.o fwi_propagator.o
	$(CC) $(DEFINES) $(CFLAGS) $(INCS) $+ -o $@ $(LIBS)

fwi_generatemodel.o:fwi_generatemodel.c
	$(CC) $(DEFINES) $(CFLAGS) $(INCS) $^ -c -o $@

fwi_main.o:fwi_main.c
	$(CC) $(DEFINES) $(CFLAGS) $(INCS) $^ -c -o $@

fwi_common.o:fwi_common.c
	$(CC) $(DEFINES) $(CFLAGS) $(INCS) $^ -c -o $@

fwi_kernel.o:fwi_kernel.c
	$(CC) $(DEFINES) $(CFLAGS) $(INCS) $^ -c -o $@

fwi_propagator.o:fwi_propagator.c
	$(CC) $(DEFINES) $(CFLAGS) $(INCS) $^ -c -o $@


.PHONY:all clean run debug memcheck input

clean:
	rm -rf *.o *.optrpt $(TARGET) fname* InputModels/* results/*

input: ModelGenerator
	./ModelGenerator ../SetupParams/fwi_params.txt ../SetupParams/fwi_frequencies.txt

run: fwi.intel64
	./fwi.intel64 ../SetupParams/fwi_params.txt ../SetupParams/fwi_frequencies.txt

debug: fwi.intel64
	gdb --args fwi.intel64 ../SetupParams/fwi_params.txt ../SetupParams/fwi_frequencies.txt

memcheck: fwi.intel64
	valgrind --leak-check=full ./fwi.intel64 ../SetupParams/fwi_params.txt ../SetupParams/fwi_frequencies.txt > valgrind.log 2>&1

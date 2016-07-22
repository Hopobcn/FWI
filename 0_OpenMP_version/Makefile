# possible values <INTEL,GNU>
COMPILER   =INTEL
# possible values <YES,NO>
PERFORM_IO =YES
# possible values <YES,NO>
DEBUG      =NO


ifeq ($(COMPILER),INTEL)
	CC        = icc
	# CFLAGS    = -std=c99 -O2 -restrict -qopenmp #-qopt-report=2
	# -opt-streaming-cache-evict=0 -opt-streaming-stores always
	CFLAGS    = -std=c99 -xHost -qopenmp -restrict -prec-div -fp-model precise -ftz -opt-prefetch-distance=64
	CINCS     =
	CLIBS     = -lm
else
	CC        = gcc
	CFLAGS    = -std=c99 -O3 -fopenmp -Wno-unused-result
	CINCS     =
	CLIBS     = -lm
endif

ifeq ($(PERFORM_IO),NO)
	DEFINES  += -DDO_NOT_PERFORM_IO
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

fwi_generatemodel.o:fwi_generatemodel.c
	$(CC) $(DEFINES) $(CFLAGS)  $(CINCS) $^ -c -o $@

fwi_main.o:fwi_main.c
	$(CC) $(DEFINES) $(CFLAGS)  $(CINCS) $^ -c -o $@

fwi_common.o:fwi_common.c
	$(CC) $(DEFINES) $(CFLAGS) $(CINCS) $^ -c -o $@

fwi_kernel.o:fwi_kernel.c
	$(CC) $(DEFINES) $(CFLAGS) $(CINCS) $^ -c -o $@

fwi_propagator.o:fwi_propagator.c
	$(CC) $(DEFINES) $(CFLAGS) $(CINCS) $^ -c -o $@


.PHONY:all clean run debug memcheck input printvars

clean:
	rm -rf *.o *.optrpt $(TARGET) *.log Results/*

input:
	./ModelGenerator ../SetupParams/fwi_params.txt ../SetupParams/fwi_frequencies.txt

run:
	./fwi.intel64 ../SetupParams/fwi_params.txt ../SetupParams/fwi_frequencies.txt

debug:
	gdb --args fwi.intel64 ../SetupParams/fwi_params.txt ../SetupParams/fwi_frequencies.txt

memcheck:
	valgrind --leak-check=full ./fwi.intel64 ../SetupParams/fwi_params.txt ../SetupParams/fwi_frequencies.txt > valgrind.log 2>&1

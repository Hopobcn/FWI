# Use -k flag to generate intermediate code
CC     = imcc
CFLAGS = --ompss --instrument --Wn,-xHost
LDFLAGS = -ldl -lrt -lpthread
#DEFINES = -DDEBUG -DDO_NOT_PERFORM_IO

#DEFINESS += -DDO_NOT_PERFORM_IO 

TARGET=fwi.intel64 ModelGenerator

all:$(TARGET)

fwi.intel64:fwi_main.o fwi_common.o fwi_kernel.o fwi_propagator.o
	$(CC) $(DEFINES) $(CFLAGS) $(INCS) $+ -o $@ $(LIBS) $(LDFLAGS)

ModelGenerator:fwi_generatemodel.o fwi_common.o fwi_kernel.o fwi_propagator.o
	$(CC) $(DEFINES) $(CFLAGS) $(INCS) $+ -o $@ $(LIBS) $(LDFLAGS)

fwi_generatemodel.o:fwi_generatemodel.c 
	$(CC) $(DEFINES) $(CFLAGS)  $(INCS) $+ -c

fwi_main.o:fwi_main.c 
	$(CC) $(DEFINES) $(CFLAGS)  $(INCS) $+ -c

fwi_common.o:fwi_common.c
	$(CC) $(DEFINES) $(CFLAGS) $(INCS) $+ -c

fwi_kernel.o:fwi_kernel.c
	$(CC) $(DEFINES) $(CFLAGS) $(INCS) $+ -c

fwi_propagator.o:fwi_propagator.c
	$(CC) $(DEFINES) $(CFLAGS) $(INCS) $+ -c


.PHONY:all clean run debug memcheck input

clean:
	rm -rf *.o *.optrpt $(TARGET) *.row *.prv *.pcf InputModels/* results/*

input:
	./ModelGenerator ../SetupParams/fwi_params.txt ../SetupParams/fwi_frequencies.txt

run:
	OMP_NUM_THREADS=16 ./fwi.intel64 ../SetupParams/fwi_params.txt ../SetupParams/fwi_frequencies.txt

irun:
	NX_ARGS="--summary" ./trace.sh ./fwi.intel64 ../SetupParams/fwi_params.txt ../SetupParams/fwi_frequencies.txt


debug:
	gdb --args fwi.intel64 ../SetupParams/fwi_params.txt ../SetupParams/fwi_frequencies.txt

memcheck:
	valgrind ./fwi.intel64 ../SetupParams/fwi_params.txt ../SetupParams/fwi_frequencies.txt

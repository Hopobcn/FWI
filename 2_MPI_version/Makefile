CC     = mpicc
CFLAGS = -std=c99 -O1 -openmp
DEFINES = -DDO_NOT_PERFOM_IO -DDEBUG
LDFLAGS = 
# podemos generar el codigo intermedio con -k

#CC = icc
#CFLAGS = -g -std=c99 -O3 -xHost -openmp -qopt-assume-safe-padding


#CFLAGS = -std=c99 -g -openmp
#CFLAGS = --opmss --Wn,-std=gnu99 --Wn,-xHost


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
	rm -rf *.o *.optrpt $(TARGET) *.row *.prv *.pcf InputModels/* results/* *.err *.out core* mpi_*

input:
	./ModelGenerator ../SetupParams/fwi_params.txt ../SetupParams/fwi_frequencies.txt

run:
	bsub < jobscript.sh

irun:
	NX_ARGS="--summary" ./trace.sh ./fwi.intel64 ../SetupParams/fwi_params.txt ../SetupParams/fwi_frequencies.txt


debug:
	gdb --args fwi.intel64 ../SetupParams/fwi_params.txt ../SetupParams/fwi_frequencies.txt

memcheck:
	valgrind ./fwi.intel64 ../SetupParams/fwi_params.txt ../SetupParams/fwi_frequencies.txt

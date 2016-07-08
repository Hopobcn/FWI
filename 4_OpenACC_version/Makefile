CC       = pgcc
CFLAGS   = -g -c99 -fast
INCS     =
LIBS     = -lm
DEFINES  += -DDO_NOT_PERFOM_IO

ifeq ($(CC),gcc)
	CFLAGS += -Wno-unused-result
endif

TARGET=fwi.intel64 ModelGenerator

all:$(TARGET)

fwi.intel64: fwi_main.o fwi_common.o fwi_kernel.o fwi_propagator.o
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

input:
	./ModelGenerator ../SetupParams/fwi_params.txt ../SetupParams/fwi_frequencies.txt

run:
	./fwi.intel64 ../SetupParams/fwi_params.txt ../SetupParams/fwi_frequencies.txt

debug:
	gdb --args fwi.intel64 ../SetupParams/fwi_params.txt ../SetupParams/fwi_frequencies.txt

memcheck:
	valgrind --leak-check=full ./fwi.intel64 ../SetupParams/fwi_params.txt ../SetupParams/fwi_frequencies.txt > valgrind.log 2>&1

# Compiler and libraries setup
MPICC    =mpimcc
CC       =icc
CFLAGS   =-Iinclude -DMPI_ENABLED -DMPICH_IGNORE_CXX_SEEK -DMPICH_SKIP_MPICXX -O3 -xHost -std=gnu99
MCFLAGS  =--ompss --no-copy-deps
LDFLAGS  =-lrt -lm -lstdc++
PRECFLAGS=--Wn,
CFLAGS   :=$(addprefix $(PRECFLAGS), $(CFLAGS) )
#DEFINES  =-DDEBUG
#DEFINES  +=-DDO_NOT_PERFOM_IO
DEFINES  +=-DUSE_NMVE

TARGET=fwi.intel64 create_schedule ModelGenerator

all:$(TARGET)

fwi.intel64:fwi_main.o fwi_common.o fwi_kernel.o fwi_propagator.o fwi_offload.o fwi_sched.o
	$(MPICC) $(DEFINES) $(CFLAGS) $(MCFLAGS) $(LDFLAGS) $^ -o $@

fwi_main.o:fwi_main.c
	$(MPICC) $(DEFINES) $(CFLAGS) $(MCFLAGS) -c $^ -o $@

fwi_common.o:fwi_common.c
	$(MPICC) $(DEFINES) $(CFLAGS) $(MCFLAGS) -c $^ -o $@

fwi_kernel.o:fwi_kernel.c
	$(MPICC) $(DEFINES) $(CFLAGS) $(MCFLAGS) -c $^ -o $@

fwi_propagator.o:fwi_propagator.c
	$(MPICC) $(DEFINES) $(CFLAGS) $(MCFLAGS) -c $^ -o $@

fwi_offload.o:fwi_offload.c
	$(MPICC) $(DEFINES) $(CFLAGS) $(MCFLAGS) -c $^ -o $@

fwi_sched.o:fwi_sched.c
	$(MPICC) $(DEFINES) $(CFLAGS) $(MCFLAGS) -c $^ -o $@

create_schedule:fwi_schedule.c
	$(CC) $(DEFINES) -std=gnu99 fwi_schedule.c -o create_schedule -lm

ModelGenerator:fwi_generatemodel.c fwi_common.o
	$(MPICC) $(LDFLAGS) -std=gnu99 fwi_generatemodel.c fwi_common.o -o ModelGenerator -lm

.PHONY:all clean input submit schedule cleanout

clean:
	rm -rf *.o *.optrpt $(TARGET) fname* InputModels/* results/*

input:
	./ModelGenerator ../SetupParams/fwi_params.txt ../SetupParams/fwi_frequencies.txt

run:
	./fwi.intel64 ../SetupParams/fwi_params.txt ../SetupParams/fwi_frequencies.txt

submit:
	qsub < sdv_x2x.pbs && watch qstat

schedule:
	./create_schedule fwi_params.txt fwi_frequencies.txt


cleanout:
	rm -rf worker.* slave.* *.err *.out core*

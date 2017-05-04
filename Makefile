# Please set MCXX_HOME & NANOS_HOME


CC=gcc
MCC=mcc
PGCC=pgcc
##################################################################################
CFLAGS=-Wall -g
MCCFLAGS=--ompss -O3
PGCFLAGS=-O3 -fast -Minline -acc -ta=tesla:cuda8.0,cc30,cc35,cc50,cc60,maxregcount:128 -Minfo=accel
##################################################################################

DEFS=-DDO_NOT_PERFORM_IO -DUSE_OMPSS

INC=-I$(NANOS_HOME)/include/nanox -Iinclude
LIBS=-L$(NANOS_HOME)/lib/performance -lnanox-ompss -lnanox-c -lnanox-gpu-api

OBJS = main/fwi_main.o src/fwi_common.o src/fwi_core.o src/fwi_kernel.o src/fwi_propagator.o

all: fwi

fwi: $(OBJS)
	$(PGCC) $(PGCFLAGS) $(DEFS) $(OBJS) -o fwi $(LIBS)

main/fwi_main.o: main/fwi_main.c
	$(MCC) $(MCCFLAGS) $(DEFS) $(INC) -o main/fwi_main.o  -c main/fwi_main.c

src/fwi_common.o: src/fwi_common.c
	$(MCC) $(MCCFLAGS) $(DEFS) $(INC) -o src/fwi_common.o -c src/fwi_common.c

src/fwi_core.o: src/fwi_core.c
	$(MCC) $(MCCFLAGS) $(DEFS) $(INC) -o src/fwi_core.o   -c src/fwi_core.c

src/fwi_kernel.o: src/fwi_kernel.c
	$(MCC) $(MCCFLAGS) $(DEFS) $(INC) -o src/fwi_kernel.o -c src/fwi_kernel.c

src/fwi_propagator.o: src/fwi_propagator.c
	$(MCC) $(MCCFLAGS) $(DEFS) $(INC) -o src/fwi_propagator.o -c src/fwi_propagator.c -kk
	$(CC)  $(CFLAGS)   $(DEFS) $(INC) -o mcc_fwi_propagator.o -c mcc_fwi_propagator.c
	$(PGCC) $(PGCFLAGS) $(DEFS) $(INC) -o openacc_fwi_propagator.o -c openacc_fwi_propagator.c
	$(LD) -r mcc_fwi_propagator.o openacc_fwi_propagator.o -o src/fwi_propagator.o

clean:
	rm -f $(OBJS) fwi *.o openacc_* mcc_* *.log

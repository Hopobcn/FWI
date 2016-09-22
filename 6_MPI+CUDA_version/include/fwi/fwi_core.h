//TODO PUT COPYRIGHT 
#ifndef _FWI_CORE_H_
#define _FWI_CORE_H_

#include "fwi_kernel.h"

void kernel( propagator_t propagator, real waveletFreq, int shotid, char* outputfolder, char* shotfolder);

void gather_shots( char* outputfolder, const real waveletFreq, const int nshots, const int numberOfCells );

int execute_simulation( int argc, char* argv[] );

#endif /* end of _FWI_CORE_H_ definition */

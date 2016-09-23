/*
 * =====================================================================================
 *
 *       Filename:  fwi_sched.c
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  04/04/2016 09:03:09 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Samuel Rodriguez Bernabeu (), samuel.rodriguez@bsc.es
 *   Organization:  Barcelona Supercomputing Center (BSC)
 *
 * =====================================================================================
 */

#ifndef _FWI_SCHED_H_
#define _FWI_SCHED_H_

#include <stdio.h>
#include <stdlib.h>
#include "fwi_common.h"


typedef struct{
	integer nfreqs;
	integer nshots;
	integer ngrads;
	integer ntests;
	char    outputfolder[200];

	real    *freq;
	integer *forws;
	integer *backs;
	integer *stacki;
	real    *dt;
	real    *dz;
	real    *dx;
	real    *dy;	
	integer    *dimmz;
	integer    *dimmx;
	integer    *dimmy;
	integer *ppd;
	integer *nworkers;
} schedule_t;

void schedule_free( schedule_t S );


schedule_t load_schedule( const char* filename ); 


#endif /*  end of _FWI_SCHED_H_ definition */

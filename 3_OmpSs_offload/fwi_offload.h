/*
 * =====================================================================================
 *
 *       Filename:  main.c
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  04/03/2016 04:58:23 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Samuel Rodriguez Bernabeu (), samuel.rodriguez@bsc.es
 *   Organization:  Barcelona Supercomputing Center (BSC)
 *
 * =====================================================================================
 */

#ifndef _FWI_OFFLOAD_H_
	#define _FWI_OFFLOAD_H_
	
	#include <stdlib.h>
	#include <stdio.h>
	#include <mpi.h>

	typedef struct{
		MPI_Comm spawn_comm;
		MPI_Comm intercomm;
		int size;
	} booster_alloc_t;

	booster_alloc_t initialize_alloc_attr( MPI_Comm comm, int hosts, int pph, int offset );

	booster_alloc_t allocate_slaves( int hosts );

	booster_alloc_t allocate_workers( int hosts, int pph, int id );


#endif /*  end of  _FWI_OFFLOAD_H_ definition */

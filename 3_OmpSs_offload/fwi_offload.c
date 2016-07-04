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
#include "fwi_offload.h"

booster_alloc_t initialize_alloc_attr( MPI_Comm comm, int hosts, int pph, int offset )
{
	booster_alloc_t attributes;
	attributes.spawn_comm = comm;
	attributes.size = hosts * pph;

	deep_booster_alloc_offset( comm, hosts, pph, &attributes.intercomm, offset );

	if ( attributes.intercomm != MPI_COMM_NULL ){
		return ( attributes );
	} else {
		fprintf(stderr, "ERROR: OmpSs offload could not allocate some nodes\n");
		abort();
	}
};

booster_alloc_t allocate_slaves( int hosts )
{
#ifdef DEBUG
	fprintf(stderr, "Master node is allocating %d slave nodes\n", hosts);
#endif

	return ( initialize_alloc_attr( MPI_COMM_WORLD, hosts, 1, 0) );
};

booster_alloc_t allocate_workers( int hosts, int pph, int id )
{
#ifdef DEBUG
	fprintf(stderr, "%d-th slave node is allocating %d worker nodes\n", id, hosts);
#endif

	return ( initialize_alloc_attr ( MPI_COMM_SELF, hosts, pph, hosts * id) );
};

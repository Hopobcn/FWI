/*
 * =====================================================================================
 *
 *       Filename:  fwi_kernel.h
 *
 *    Description:  Kernel propagator implementation
 *
 *        Version:  1.0
 *        Created:  14/12/15 12:10:24
 *       Revision:  none
 *       Compiler:  icc
 *
 *         Author:  YOUR NAME (),
 *   Organization:
 *
 * =====================================================================================
 */

#ifndef _FWI_KERNEL_H_
#define _FWI_KERNEL_H_

#include "fwi_propagator.h"

/* -------------------- MEMORY RELATED FUNCTIONS ------------------------------ */

void alloc_memory_shot( const integer numberOfCells,
                        coeff_t *c,
                        s_t     *s,
                        v_t     *v,
											  real    **rho);

void free_memory_shot( coeff_t *c,
                       s_t     *s,
                       v_t     *v,
											 real    **rho);

void check_memory_shot( const integer numberOfCells,
                        coeff_t *c,
                        s_t     *s,
                        v_t     *v,
												real    *rho);

/* --------------- I/O RELATED FUNCTIONS -------------------------------------- */

void load_initial_model ( const real    waveletFreq,
													const integer dimmz,
													const integer dimmx,
													const integer dimmy,
                         	coeff_t *c,
                         	s_t     *s,
                         	v_t     *v,
                         	real    *rho);

void write_snapshot ( char          *folder,
                      const int     suffix,
                      real          *data,
                      const integer numberOfCells,
                      const int     nfields);

void read_snapshot ( char          *folder,
                     const int     suffix,
                     real          *data,
                     const integer numberOfCells,
                     const int     nfields);


/* --------------- WAVE PROPAGATOR FUNCTIONS --------------------------------- */

void propagate_shot ( time_d        direction,
                     v_t           v,
                     s_t           s,
                     coeff_t       coeffs,
                     real          *rho,
                     int           timesteps,
                     int           ntbwd,
                     real          dt,
                     real          dzi,
                     real          dxi,
                     real          dyi,
                     integer       nz0,
                     integer       nzf,
                     integer       nx0,
                     integer       nxf,
                     integer       ny0,
                     integer       nyf,
                     integer       stacki,
                     char          *folder,
                     real          *dataflush,
                     integer       datalen,
                     integer       dimmz,
                     integer       dimmx);


/* --------------- BOUNDARY EXCHANGES ---------------------------------------- */

/*
NAME:exchange_boundaries
PURPOSE: data exchanges between the boundary layers of the analyzed volume

v                   (in) struct containing velocity arrays (4 points / cell x 3 components / point = 12 arrays)
numElement          (in) Number of elements to exchange
rank                (in) MPI process identifier
numTasks            (in) MPI number of tasks
idxt                (in) identifier related to the folder
nyf                 (in) final plane to be exchanged
ny0                 (in) intial plane to be exchanged

RETURN none
*/

inline integer exchange_buffer ( real* _bufferA, real* _bufferB, integer _torank, uint64_t message_size)
{
	MPI_Status stat;
	int tag = 100;

	MPI_Sendrecv( &_bufferA, message_size, MPI_FLOAT, _torank, tag,
                &_bufferB, message_size, MPI_FLOAT, _torank, tag,
               	MPI_COMM_WORLD, &stat);
};

void exchange_velocity_boundaries ( v_t *v, 
																		integer numberOfCells, 
																		integer nyf, 
																		integer ny0 );

void exchange_stress_boundaries ( s_t *s, 
																	integer numberOfCells, 
																	integer nyf, 
																	integer ny0 );


#endif /* end of _FWI_KERNEL_H_ definition */

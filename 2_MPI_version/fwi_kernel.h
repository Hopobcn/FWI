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

#define EXCHANGE(sendbuf, recvbuf, dst, src, count) {                             \
    exchange_buffer((sendbuf),(recvbuf),(dst),(src),(count), __FILE__, __LINE__); \
}                                                                                 
inline integer exchange_buffer (const real*   sendbuf, 
                                      real*   recvbuf, 
                                const integer dst, 
                                const integer src, 
                                const integer message_size,
                                const char*   file,
                                const integer line)
{
    int tag = 100;
#if defined(DEBUG) && defined(DEBUG_MPI)
    log_info( "         [BEFORE]MPI sendrecv [count:%d][dst:%d][src:%d] %s : %d", message_size,  dst, src, file, line);
#endif

#if 0
    MPI_Status status;
    //MPI_Sendrecv may deadlock in some MPI implementations!!!!
    //            --> (1) order communications or (2) use non-blocking calls
    MPI_Sendrecv( sendbuf, message_size, MPI_FLOAT, dst, tag,
                  recvbuf, message_size, MPI_FLOAT, src, tag,
                  MPI_COMM_WORLD, &status);
#else
    MPI_Status  statuses[2];
    MPI_Request requests[2];
    
    MPI_Irecv( recvbuf, message_size, MPI_FLOAT, dst, tag, MPI_COMM_WORLD, &requests[0] );
    MPI_Isend( sendbuf, message_size, MPI_FLOAT, dst, tag, MPI_COMM_WORLD, &requests[1] );
    MPI_Waitall(2, requests, statuses);
#endif

#if defined(DEBUG) && defined(DEBUG_MPI)
    log_info( "         [AFTER ]MPI sendrecv                          %s : %d", file, line);
#endif
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

/*
 * =============================================================================
 * Copyright (c) 2016-2018, Barcelona Supercomputing Center (BSC)
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *     * Redistributions of source code must retain the above copyright
 *       notice, this list of conditions and the following disclaimer.
 *     * Redistributions in binary form must reproduce the above copyright
 *       notice, this list of conditions and the following disclaimer in the
 *       documentation and/or other materials provided with the distribution.
 *     * Neither the name of the <organization> nor the
 *       names of its contributors may be used to endorse or promote products
 *       derived from this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL <COPYRIGHT HOLDER> BE LIABLE FOR ANY
 * DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 * =============================================================================
 */

#ifndef _FWI_KERNEL_H_
#define _FWI_KERNEL_H_

#include "fwi_propagator.h"

/*
 * Ensures that the domain contains a minimum number of planes.
 * Just needed for debugging when running small cases.
 */
void check_domain_dimensions ( const integer dimmz,
                               const integer dimmx,
                               const integer dimmy);

void set_array_to_random_real(real* restrict array,
                              const integer length);

void set_array_to_constant(real* restrict array,
                           const real value,
                           const integer length);

void alloc_memory_shot( const integer dimmz,
                        const integer dimmx,
                        const integer dimmy,
                        coeff_t *c,
                        s_t     *s,
                        v_t     *v,
                        real    **rho);

void free_memory_shot( coeff_t *c,
                       s_t     *s,
                       v_t     *v,
                       real    **rho);

void check_memory_shot( const integer dimmz,
                        const integer dimmx,
                        const integer dimmy,
                        coeff_t *c,
                        s_t     *s,
                        v_t     *v,
                        real    *rho);

/* --------------- I/O RELATED FUNCTIONS -------------------------------------- */

void load_local_velocity_model ( const real    waveletFreq,
                                 const integer dimmz,
                                 const integer dimmx,
                                 const integer FirstYPlane,
                                 const integer LastYPlane,
                                 coeff_t *c,
                                 s_t     *s,
                                 v_t     *v,
                                 real    *rho);

void write_snapshot ( char         *folder,
                      const int     suffix,
                      v_t          *v,
                      const integer dimmz,
                      const integer dimmx,
                      const integer dimmy);

void read_snapshot ( char         *folder,
                     const int     suffix,
                     v_t          *v,
                     const integer dimmz,
                     const integer dimmx,
                     const integer dimmy);


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
                     real          *UNUSED(dataflush),
                     integer       dimmz,
                     integer       dimmx,
                     integer       dimmy);


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

#if defined(USE_MPI)
#define EXCHANGE(sendbuf, recvbuf, dst, src, count) {                             \
    exchange_buffer((sendbuf),(recvbuf),(dst),(src),(count), __FILE__, __LINE__); \
}                                                                                 
static inline integer exchange_buffer (const real*   sendbuf, 
                                             real*   recvbuf, 
                                       const integer dst, 
                                       const integer src, 
                                       const integer message_size,
                                       const char*   file,
                                       const integer line)
{
    int err;
    int tag = 100;
    
    print_debug( "         [BEFORE]MPI sendrecv [count:%d][dst:%d][src:%d] %s : %d", message_size,  dst, src, file, line);

    MPI_Status  statuses[2];
    MPI_Request requests[2];

#if defined(_OPENACC)
    #pragma acc host_data use_device(recvbuf, sendbuf)
#endif
    {
        MPI_Irecv( recvbuf, message_size, MPI_FLOAT, dst, tag, MPI_COMM_WORLD, &requests[0] );
        MPI_Isend( sendbuf, message_size, MPI_FLOAT, dst, tag, MPI_COMM_WORLD, &requests[1] );
    }

    err = MPI_Waitall(2, requests, statuses);

    print_debug( "         [AFTER ]MPI sendrecv                          %s : %d", file, line);
    
    return err;
};

void exchange_velocity_boundaries ( v_t v, 
                                    const integer plane_size, 
                                    const integer nyf, 
                                    const integer ny0 );

void exchange_stress_boundaries ( s_t s, 
                                  const integer plane_size, 
                                  const integer nyf, 
                                  const integer ny0 );
#else
#define EXCHANGE(sendbuf, recvbuf, dst, src, count)
#endif


#endif /* end of _FWI_KERNEL_H_ definition */

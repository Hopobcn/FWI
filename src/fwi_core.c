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

#include "fwi/fwi_core.h"
#include "fwi/fwi_sched.h"

/*
 * In order to generate a source for injection,
 * /system/support/bscgeo/src/wavelet.c
 * functions can be used.
 */
void kernel( propagator_t propagator, real waveletFreq, int shotid, char* outputfolder, char* shotfolder)
{
#if defined(USE_MPI)
    /* find ourselves into the MPI space */
    int mpi_rank, mpi_size;
    MPI_Comm_rank( MPI_COMM_WORLD, &mpi_rank);
    MPI_Comm_size( MPI_COMM_WORLD, &mpi_size);
#endif /* USE_MPI */

    /* local variables */
    int stacki;
    double start_t, end_t;
    real dt,dz,dx,dy;
    integer dimmz, dimmx, dimmy, MaxYPlanesPerWorker, forw_steps, back_steps;

    load_shot_parameters( shotid, &stacki, &dt, &forw_steps, &back_steps,
            &dz, &dx, &dy,
            &dimmz, &dimmx, &dimmy,
            &MaxYPlanesPerWorker,
            outputfolder, waveletFreq );

#if defined(USE_MPI)
    /* aux variables, just to make it more readable */
    const int FIRSTRANK = 0;
    const int LASTRANK  = mpi_size - 1;

    /* Compute the integration limits in order to load the correct slice from the input
     * velocity model. These are not the limits for the wave propagator! (they are local,
     * i.e. starts at zero!) */
    const integer y0 = (mpi_rank == FIRSTRANK) ? 0     : (MaxYPlanesPerWorker * mpi_rank) - HALO;
    const integer yf = (mpi_rank == LASTRANK ) ? dimmy : y0 + MaxYPlanesPerWorker;
    const integer edimmy = (yf - y0);
#else
    const integer y0 = 0;
    const integer yf = dimmy;
    const integer edimmy = dimmy;
#endif /* USE_MPI */

    /* Compute integration limits for the wave propagator. 
     * It assumes that the volume is local, so the indices start at zero */
    const integer nz0 = 0;
    const integer ny0 = 0;
    const integer nx0 = 0;
    const integer nzf = dimmz;
    const integer nxf = dimmx;
    const integer nyf = edimmy;
    const integer numberOfCells = dimmz * dimmx * edimmy;

    real    *rho;
    v_t     v;
    s_t     s;
    coeff_t coeffs;

    print_debug("The length of local arrays is " I " cells zxy[%d][%d][%d]", numberOfCells, nzf, nxf, nyf);

    /* allocate shot memory */
    alloc_memory_shot  ( dimmz, dimmx, (nyf - ny0), &coeffs, &s, &v, &rho);

    /* load initial model from a binary file */
    load_local_velocity_model ( waveletFreq, dimmz, dimmx, y0, yf, &coeffs, &s, &v, rho);

    /* Allocate memory for IO buffer */
    real* io_buffer = (real*) __malloc( ALIGN_REAL, numberOfCells * sizeof(real) * WRITTEN_FIELDS );

    /* inspects every array positions for leaks. Enabled when DEBUG flag is defined */
    check_memory_shot  ( dimmz, dimmx, (nyf - ny0), &coeffs, &s, &v, rho);

    /* Perform forward, backward or test propagations */
    switch( propagator )
    {
    case( RTM_KERNEL ):
    {
        start_t = dtime();

        propagate_shot ( FORWARD,
                         v, s, coeffs, rho,
                         forw_steps, back_steps -1,
                         dt,dz,dx,dy,
                         nz0, nzf, nx0, nxf, ny0, nyf,
                         stacki,
                         shotfolder,
                         io_buffer,
                         dimmz, dimmx, (nyf - ny0));

        end_t = dtime();

        print_stats("Forward propagation finished in %lf seconds", end_t - start_t );

        start_t = dtime();

        propagate_shot ( BACKWARD,
                         v, s, coeffs, rho,
                         forw_steps, back_steps -1,
                         dt,dz,dx,dy,
                         nz0, nzf, nx0, nxf, ny0, nyf,
                         stacki,
                         shotfolder,
                         io_buffer,
                         dimmz, dimmx, (nyf - ny0));

        end_t = dtime();

        print_stats("Backward propagation finished in %lf seconds", end_t - start_t );

#if defined(DO_NOT_PERFORM_IO)
        print_info("Warning: we are not creating gradient nor preconditioner "
                   "fields, because IO is not enabled for this execution" );
#else

#if defined(USE_MPI)
        if ( mpi_rank == 0 ) 
#endif /* USE_MPI */
        {
            char fnameGradient[300];
            char fnamePrecond[300];
            sprintf( fnameGradient, "%s/gradient_%05d.dat", shotfolder, shotid );
            sprintf( fnamePrecond , "%s/precond_%05d.dat" , shotfolder, shotid );

            FILE* fgradient = safe_fopen( fnameGradient, "wb", __FILE__, __LINE__ );
            FILE* fprecond  = safe_fopen( fnamePrecond , "wb", __FILE__, __LINE__ );

            print_info("Storing local preconditioner field in %s", fnameGradient );
            safe_fwrite( io_buffer, sizeof(real), numberOfCells * 12, fgradient, __FILE__, __LINE__ );

            print_info("Storing local gradient field in %s", fnamePrecond);
            safe_fwrite( io_buffer, sizeof(real), numberOfCells * 12, fprecond , __FILE__, __LINE__ );

            safe_fclose( fnameGradient, fgradient, __FILE__, __LINE__ );
            safe_fclose( fnamePrecond , fprecond , __FILE__, __LINE__ );
        }
#endif /* end DO_NOT_PERFORM_IO */

        break;
    }
    case( FM_KERNEL  ):
    {
        start_t = dtime();

        propagate_shot ( FWMODEL,
                         v, s, coeffs, rho,
                         forw_steps, back_steps -1,
                         dt,dz,dx,dy,
                         nz0, nzf, nx0, nxf, ny0, nyf,
                         stacki,
                         shotfolder,
                         io_buffer,
                         dimmz, dimmx, dimmy);

        end_t = dtime();

        print_stats("Forward Modelling finished in %lf seconds", end_t - start_t );
       
        break;
    }
    default:
    {
        print_error("Invalid propagation identifier");
        abort();
    }
    } /* end case */

    // liberamos la memoria alocatada en el shot
    free_memory_shot  ( &coeffs, &s, &v, &rho);
    __free( io_buffer );
};

void gather_shots( char* outputfolder, const real waveletFreq, const int nshots, const int numberOfCells )
{
#if defined(DO_NOT_PERFORM_IO)
    print_info("Warning: we are not gathering the results because the IO is disabled "
               "for this execution");
#else
    /* ---------  GLOBAL PRECONDITIONER ACCUMULATION --------- */
    print_info("Gathering local preconditioner fields");

    /* variables for timming */
    double start_t, end_t;

    /* buffers to read and accumulate the fields */
    real* sumbuffer  = (real*)  __malloc( ALIGN_REAL, numberOfCells * sizeof(real) * WRITTEN_FIELDS ); 
    real* readbuffer = (real*)  __malloc( ALIGN_REAL, numberOfCells * sizeof(real) * WRITTEN_FIELDS );
    
    start_t = dtime();

    /* set buffer positions to zero */
    memset ( sumbuffer, 0, numberOfCells * sizeof(real) * WRITTEN_FIELDS );

    for( int shot=0; shot < nshots; shot++)
    {
        char readfilename[300];
        sprintf( readfilename, "%s/shot.%2.1f.%05d/precond_%05d.dat", 
                outputfolder, waveletFreq, shot, shot);

        print_info("Reading preconditioner file '%s'", readfilename );

        FILE* freadfile = safe_fopen( readfilename, "rb", __FILE__, __LINE__ );
        safe_fread ( readbuffer, sizeof(real), numberOfCells * WRITTEN_FIELDS, freadfile, __FILE__, __LINE__ );

#if defined(_OPENMP)
        #pragma omp parallel for
#endif
#if defined(__INTEL_COMPILER)
        #pragma simd
#endif
        for( int i = 0; i < numberOfCells * WRITTEN_FIELDS; i++)
            sumbuffer[i] += readbuffer[i];

        fclose (freadfile);
    }

    char precondfilename[300];
    sprintf( precondfilename, "%s/Preconditioner.%2.1f", outputfolder, waveletFreq );
    FILE* precondfile = safe_fopen( precondfilename, "wb", __FILE__, __LINE__ );
    safe_fwrite ( sumbuffer, sizeof(real), numberOfCells * WRITTEN_FIELDS, precondfile, __FILE__, __LINE__ );
    safe_fclose( precondfilename, precondfile, __FILE__, __LINE__ );

    end_t = dtime();

    print_stats("Gatering process for preconditioner %s (freq %2.1f) " 
                "completed in: %lf seconds",  
                precondfilename, waveletFreq, end_t - start_t  );

    /* ---------  GLOBAL GRADIENT ACCUMULATION --------- */
    print_info("Gathering local gradient fields");

    start_t = dtime();

    /* set buffer positions to zero */
    memset ( sumbuffer, 0, numberOfCells * sizeof(real) * WRITTEN_FIELDS );

    for( int shot=0; shot < nshots; shot++)
    {
        char readfilename[300];
        sprintf( readfilename, "%s/shot.%2.1f.%05d/gradient_%05d.dat", 
                outputfolder, waveletFreq, shot, shot);

        print_info("Reading gradient file %s", readfilename );

        FILE* freadfile = safe_fopen( readfilename, "rb", __FILE__, __LINE__ );
        safe_fread ( readbuffer, sizeof(real), numberOfCells * WRITTEN_FIELDS, freadfile, __FILE__, __LINE__ );

#if defined(_OPENMP)
        #pragma omp parallel for
#endif
#ifdef __INTEL_COMPILER
        #pragma simd
#endif
        for( int i = 0; i < numberOfCells * WRITTEN_FIELDS; i++)
            sumbuffer[i] += readbuffer[i];

        fclose (freadfile);
    }

    char gradientfilename[300];
    sprintf( gradientfilename, "%s/Gradient.%2.1f", outputfolder, waveletFreq );
    FILE* gradientfile = safe_fopen( gradientfilename, "wb", __FILE__, __LINE__ );
    safe_fwrite ( sumbuffer, sizeof(real), numberOfCells * WRITTEN_FIELDS, gradientfile, __FILE__, __LINE__ );
    safe_fclose( gradientfilename, gradientfile, __FILE__, __LINE__ );

    end_t = dtime();

    print_stats("Gatering process for gradient %s (freq %2.1f) "        
                "completed in: %lf seconds", 
                precondfilename, waveletFreq, end_t - start_t  );

    __free(  sumbuffer);
    __free( readbuffer);
#endif /* end DO_NOT_PERFORM_IO */
};

int execute_simulation( int argc, char* argv[] )
{
#if defined(USE_MPI)
    MPI_Init ( &argc, &argv );
    int mpi_rank;
    MPI_Comm_rank( MPI_COMM_WORLD, &mpi_rank);
#elif !defined(USE_MPI) && defined(_OPENACC)
    //TODO: fix name
    int mpi_rank = 0;
#endif

#if defined(_OPENACC)
    acc_init(acc_device_default);
    int gpuid = mpi_rank % acc_get_num_devices( acc_device_default );
    acc_set_device_num( gpuid, acc_device_default );
    fprintf(stdout, "MPI rank %d with GPU %d (%d)\n", 
            mpi_rank, acc_get_device_num(acc_device_default), acc_get_num_devices(acc_device_default));
#endif /*_OPENACC*/


    /* Load parameters from schedule file */
    schedule_t s = load_schedule(argv[1]);

    for(int i=0; i<s.nfreqs; i++)
    {
        /* Process one frequency at a time */
        real waveletFreq   = s.freq[i];
        integer stacki     = s.stacki[i];
        real dt            = s.dt[i];
        integer forw_steps = s.forws[i];
        integer back_steps = s.backs[i];
        real dx            = s.dx[i];
        real dy            = s.dy[i];
        real dz            = s.dz[i];
        integer dimmz      = s.dimmz[i];
        integer dimmx      = s.dimmx[i];
        integer dimmy      = s.dimmy[i];
        //integer nworkers = s.nworkers[i];
        integer MaxYPlanesPerWorker = s.ppd[i];

        print_info("\n------ Computing %d-th frequency (%.2fHz). ------\n", i, waveletFreq);

        const integer numberOfCells = dimmz * dimmx * dimmx;
        const size_t VolumeMemory  = numberOfCells * sizeof(real) * 58;

        print_stats("Local domain size for freq %f [%d][%d][%d] is %lu bytes (%lf GB)", 
                    waveletFreq, dimmz, dimmx, dimmy, VolumeMemory, TOGB(VolumeMemory) );

        for(int grad=0; grad<s.ngrads; grad++) /* backward iteration */
        {
            print_info("Processing %d-gradient iteration", grad);

            for(int shot=0; shot<s.nshots; shot++)
            {
                char shotfolder[512];
                sprintf(shotfolder, "%s/shot.%2.2fHz.%03d", s.outputfolder, waveletFreq, shot);

#if defined(USE_MPI)
                if ( mpi_rank == 0 )
#endif
                {
                    create_folder( shotfolder );

                    store_shot_parameters( shot, &stacki, &dt, &forw_steps, &back_steps,
                                           &dz, &dx, &dy,
                                           &dimmz, &dimmx, &dimmy,
                                           &MaxYPlanesPerWorker,
                                           s.outputfolder, waveletFreq );
                }
#if defined(USE_MPI)
                MPI_Barrier( MPI_COMM_WORLD );
#endif

                kernel( RTM_KERNEL, waveletFreq, shot, s.outputfolder, shotfolder);

                print_info("\tGradient loop processed for %d-th shot", shot);

                //update_shot()
            }

//#if defined(USE_MPI)
//           MPI_Barrier( MPI_COMM_WORLD );
//
//           if ( mpi_rank == 0 ) { 
//               gather_shots( outputfolder, waveletFreq, nshots, numberOfCells );    
//           }
//
//           MPI_Barrier( MPI_COMM_WORLD );
//#else
//           gather_shots( s.outputfolder, waveletFreq, s.nshots, numberOfCells );
//#endif

            for(int test=0; test<s.ntests; test++)
            {
                print_info("\tProcessing %d-th test iteration", test);

                for(int shot=0; shot<s.nshots; shot++)
                {
                    char shotfolder[512];
                    sprintf(shotfolder, "%s/test.%05d.shot.%2.2fHz.%03d", 
                            s.outputfolder, test, waveletFreq, shot);

#if defined(USE_MPI)
                    if ( mpi_rank == 0)
#endif
                    {
                        create_folder( shotfolder );

                        store_shot_parameters( shot, &stacki, &dt, &forw_steps, &back_steps,
                                               &dz, &dx, &dy,
                                               &dimmz, &dimmx, &dimmy,
                                               &MaxYPlanesPerWorker,
                                               s.outputfolder, waveletFreq );
                    }
#if defined(USE_MPI)
                    MPI_Barrier( MPI_COMM_WORLD );
#endif
                    kernel( FM_KERNEL , waveletFreq, shot, s.outputfolder, shotfolder);

                    print_info("\t\tTest loop processed for the %d-th shot", shot);
                }
            } /* end of test loop */
        } /* end of gradient loop */
    } /* end of frequency loop */

#if defined(USE_MPI)
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();
#endif

    return 0;
}

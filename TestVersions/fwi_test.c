/*
 * =====================================================================================
 *
 *       Filename:  fwi_main.c
 *
 *    Description:  Main file of the FWI mockup
 *
 *        Version:  1.0
 *        Created:  10/12/15 10:33:40
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  YOUR NAME (),
 *   Organization:
 *
 * =====================================================================================
 */

#include "fwi_kernel.h"

/*
 * In order to generate a source for injection,
 * /system/support/bscgeo/src/wavelet.c
 * functions can be used.
 */

void kernel( propagator_t propagator, real waveletFreq, int shotid, 
             char* ref_outfolder, char* opt_outfolder, 
             char* ref_shotfolder, char* opt_shotfolder )
{
    int stacki_ref;
    real dt_ref, dz_ref, dx_ref, dy_ref;
    integer dimmz_ref, dimmx_ref, dimmy_ref;
    int forw_steps_ref, back_steps_ref;

    int stacki_opt;
    real dt_opt, dz_opt, dx_opt, dy_opt;
    integer dimmz_opt, dimmx_opt, dimmy_opt;
    int forw_steps_opt, back_steps_opt;

    load_shot_parameters( shotid, &stacki_ref, &dt_ref, &forw_steps_ref, &back_steps_ref, &dz_ref, &dx_ref, &dy_ref, &dimmz_ref, &dimmx_ref, &dimmy_ref, ref_outfolder, waveletFreq );
    load_shot_parameters( shotid, &stacki_opt, &dt_opt, &forw_steps_opt, &back_steps_opt, &dz_opt, &dx_opt, &dy_opt, &dimmz_opt, &dimmx_opt, &dimmy_opt, opt_outfolder, waveletFreq );

    // Compare reference parameters with optimized version
    fprintf(stderr, "-----: stacki\tdt\tdz\tdx\tdy\t     dimmz\tdimmx\tdimmy\tfw_steps bw_steps\n");
    fprintf(stderr, "--ref: %d     %f %f %f %f\t%d\t%d\t%d\t%d\t%d\n",
            stacki_ref, dt_ref, dz_ref, dx_ref, dy_ref, dimmz_ref, dimmx_ref, dimmy_ref, forw_steps_ref, back_steps_ref);
    fprintf(stderr, "--opt: %d     %f %f %f %f\t%d\t%d\t%d\t%d\t%d\n",
            stacki_opt, dt_opt, dz_opt, dx_opt, dy_opt, dimmz_opt, dimmx_opt, dimmy_opt, forw_steps_opt, back_steps_opt);

    assert(stacki_ref == stacki_opt);
    assert(dimmz_ref  == dimmz_opt );
    assert(dimmx_ref  == dimmz_ref );
    assert(dimmy_ref  == dimmy_ref );
    assert(forw_steps_ref == forw_steps_opt);
    assert(back_steps_ref == back_steps_opt);
    assert_eq_scalar(dz_ref, dz_opt, "dz_ref", "dz_opt", 1);
    assert_eq_scalar(dx_ref, dx_opt, "dx_ref", "dx_opt", 1);
    assert_eq_scalar(dy_ref, dy_opt, "dy_ref", "dy_opt", 1);
    assert_eq_scalar(dt_ref, dt_opt, "dt_ref", "dt_opt", 1);

    const integer numberOfCells = dimmz_ref * dimmx_ref * dimmy_ref;

    /* set LOCAL integration limits */
    const integer nz0 = 0;
    const integer ny0 = 0;
    const integer nx0 = 0;
    const integer nzf = dimmz_ref;
    const integer nxf = dimmx_ref;
    const integer nyf = dimmy_ref;

    real    *rho_ref, *rho_opt;
    v_t     v_ref, v_opt;
    s_t     s_ref, s_opt;
    coeff_t coeffs_ref, coeffs_opt;

    /* allocate shot memory */
    alloc_memory_shot  ( numberOfCells, &coeffs_ref, &s_ref, &v_ref, &rho_ref);
    alloc_memory_shot  ( numberOfCells, &coeffs_opt, &s_opt, &v_opt, &rho_opt);

    /* load initial model from a binary file */
    //load_initial_model ( waveletFreq, numberOfCells, &coeffs, &s, &v, rho);

    /* Allocate memory for IO buffer */
    //real* io_buffer = (real*) __malloc( ALIGN_REAL, numberOfCells * sizeof(real) * WRITTEN_FIELDS );

    /* inspects every array positions for leaks. Enabled when DEBUG flag is defined */
    //check_memory_shot  ( numberOfCells, &coeffs, &s, &v, rho);

    /* some variables for timming */
    double start_t, end_t;

    switch( propagator )
    {
    case( RTM_KERNEL ):
    {
        start_t = dtime();

        propagate_shot ( FORWARD,
                         v_ref, v_opt, 
                         s_ref, coeffs_ref, rho_ref,
                         forw_steps_ref, back_steps_ref -1,
                         dt_ref, dz_ref, dx_ref, dy_ref,
                         nz0, nzf, nx0, nxf, ny0, nyf,
                         stacki_ref,
                         ref_shotfolder, opt_shotfolder,
                         numberOfCells,
                         dimmz_ref, dimmx_ref);

        end_t = dtime();

        fprintf(stdout, "Forward propagation finished in %lf seconds\n", \
                         end_t - start_t );
#if 0
        start_t = dtime();
        propagate_shot ( BACKWARD,
                         v, s, coeffs, rho,
                         forw_steps, back_steps -1,
                         dt,dz,dx,dy,
                         nz0, nzf, nx0, nxf, ny0, nyf,
                         stacki,
                         shotfolder,
                         io_buffer,
                         numberOfCells,
                         dimmz, dimmx);

        end_t = dtime();

        fprintf(stdout, "Backward propagation finished in %lf seconds\n", \
                         end_t - start_t );
#endif

#if 0
        //TODO: COMPARE Gradient & Precond of ref vs opt
        char ref_fnameGradient[300], opt_fnameGradient[300];
        char ref_fnamePrecond[300],  opt_fnamePrecond[300];
        sprintf( fnameGradient, "%s/gradient_%05d.dat", shotfolder, shotid );
        sprintf( fnamePrecond , "%s/precond_%05d.dat" , shotfolder, shotid );

        FILE* fgradient = safe_fopen( fnameGradient, "wb", __FILE__, __LINE__ );
        FILE* fprecond  = safe_fopen( fnamePrecond , "wb", __FILE__, __LINE__ );

        fprintf(stderr, "Storing local preconditioner field in %s\n", fnameGradient );
        safe_fwrite( io_buffer, sizeof(real), numberOfCells * 12, fgradient, __FILE__, __LINE__ );

        fprintf(stderr, "Storing local gradient field in %s\n", fnamePrecond);
        safe_fwrite( io_buffer, sizeof(real), numberOfCells * 12, fprecond , __FILE__, __LINE__ );

        safe_fclose( fnameGradient, fgradient, __FILE__, __LINE__ );
        safe_fclose( fnamePrecond , fprecond , __FILE__, __LINE__ );
#endif
        break;
    }
    case( FM_KERNEL  ):
    {
        start_t = dtime();

        propagate_shot ( FWMODEL,
                         v_ref, v_opt,
                         s_ref, coeffs_ref, rho_ref,
                         forw_steps_ref, back_steps_ref -1,
                         dt_ref, dz_ref, dx_ref, dy_ref,
                         nz0, nzf, nx0, nxf, ny0, nyf,
                         stacki_ref,
                         ref_shotfolder, opt_shotfolder,
                         numberOfCells,
                         dimmz_ref, dimmx_ref);

        end_t = dtime();

        fprintf(stdout, "Forward Modelling finished in %lf seconds\n",  \
                         end_t - start_t );
       
        break;
    }
    default:
    {
        fprintf(stderr, "Invalid propagation identifier\n");
        abort();
    }
    } /* end case */

    // liberamos la memoria alocatada en el shot
    free_memory_shot  ( &coeffs_ref, &s_ref, &v_ref, &rho_ref);
    free_memory_shot  ( &coeffs_opt, &s_opt, &v_opt, &rho_opt);

    //__free( io_buffer );

    fprintf(stderr, "Shot memory free'd\n");
};

void gather_shots( char* outputfolder, const real waveletFreq, const int nshots, const int numberOfCells )
{
#ifdef DO_NOT_PERFORM_IO
    fprintf(stderr, "Warning: we are not doing any IO here (%s)\n", __FUNCTION__ );
#else
    /* ---------  GLOBAL PRECONDITIONER ACCUMULATION --------- */
    fprintf(stderr, "Gathering local preconditioner fields\n");

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

        fprintf(stderr, "Reading preconditioner file %s\n", readfilename );

        FILE* freadfile = safe_fopen( readfilename, "rb", __FILE__, __LINE__ );
        safe_fread ( readbuffer, sizeof(real), numberOfCells * WRITTEN_FIELDS, freadfile, __FILE__, __LINE__ );

        #pragma omp parallel for
#ifdef __INTEL_COMPILER
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

    fprintf(stderr, "Gatering process for preconditioner %s (freq %2.1f) "
                    "completed in: %lf seconds\n",
                    precondfilename, waveletFreq, end_t - start_t  );



    /* ---------  GLOBAL GRADIENT ACCUMULATION --------- */
    fprintf(stderr, "Gathering local gradient fields\n");

    start_t = dtime();

    /* set buffer positions to zero */
    memset ( sumbuffer, 0, numberOfCells * sizeof(real) * WRITTEN_FIELDS );

    for( int shot=0; shot < nshots; shot++)
    {
        char readfilename[300];
        sprintf( readfilename, "%s/shot.%2.1f.%05d/gradient_%05d.dat", 
                outputfolder, waveletFreq, shot, shot);

        fprintf(stderr, "Reading gradient file %s\n", readfilename );

        FILE* freadfile = safe_fopen( readfilename, "rb", __FILE__, __LINE__ );
        safe_fread ( readbuffer, sizeof(real), numberOfCells * WRITTEN_FIELDS, freadfile, __FILE__, __LINE__ );

        #pragma omp parallel for
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

    fprintf(stderr, "Gatering process for gradient %s (freq %2.1f) " 
                    "completed in: %lf seconds\n",    
                    precondfilename, waveletFreq, end_t - start_t );

    __free(  sumbuffer);
    __free( readbuffer);
#endif
};

int main(int argc, const char* argv[])
{
    real lenz,lenx,leny,vmin,srclen,rcvlen;
    char outputfolder[200];

    if (argc < 5) {
        fprintf(stderr, "Usage: %s [FWI_params] [FWI_freqs] [Ref.Ver. Path] [Opti.Ver. Path]\n", argv[0]);
        exit(-1);
    }

    read_fwi_parameters( argv[1], &lenz, &lenx, &leny, &vmin, &srclen, &rcvlen, outputfolder);

    const int nshots = 2;
    const int ngrads = 1;
    const int ntest  = 1;

    int   nfreqs;
    real *frequencies;

    load_freqlist( argv[2], &nfreqs, &frequencies );

    for(int i=0; i<nfreqs; i++)
    {
        /* Process one frequency at a time */
        real waveletFreq = frequencies[i];

        /* Deltas of space, 16 grid point per Hz */
        real dx = vmin / (16.0 * waveletFreq);
        real dy = vmin / (16.0 * waveletFreq);
        real dz = vmin / (16.0 * waveletFreq);

        /* number of cells along axis, adding HALO planes */
        integer dimmz = roundup(ceil( lenz / dz ) + 2*HALO, HALO);
        integer dimmy = roundup(ceil( leny / dy ) + 2*HALO, HALO);
        integer dimmx = roundup(ceil( lenx / dx ) + 2*HALO, HALO);

        fprintf(stderr, "Domain dimensions (dimm x, y, z) %d %d %d delta of space (dx,dy,dz) %f %f %f vim %f\n",
                         dimmx, dimmy, dimmz, dx, dy, dz, vmin); 

        /* compute delta T */
        real dt = 68e-6 * dx;

        /* dynamic I/O */
        int stacki = floor(  0.25 / (2.5 * waveletFreq * dt) );

        fprintf(stderr, "Stack(i) vale is %d\n", stacki);

        const integer numberOfCells = dimmz * dimmx * dimmx;
        const integer VolumeMemory  = numberOfCells * sizeof(real) * 58;

        fprintf(stderr, "Local domain size is " I " bytes (%f GB)\n", \
                         VolumeMemory, TOGB(VolumeMemory) );

        /* compute time steps */
        int forw_steps = max_int ( IT_FACTOR * (srclen/dt), 1);
        int back_steps = max_int ( IT_FACTOR * (rcvlen/dt), 1);

        fprintf(stderr, "stacki value is %d. "              \
                        "forward propagation steps: %d "    \
                        "backward propagation steps: %d\n", \
                        stacki, forw_steps, back_steps);

        for(int grad=0; grad<ngrads; grad++) /* iteracion de inversion */
        {
            fprintf(stderr, "Processing %d-th gradient iteration.\n", grad);

            for(int shot=0; shot<nshots; shot++)
            {
                char ref_shotfolder[200], opt_shotfolder[200];
                char ref_outfolder[200],  opt_outfolder[200];
                sprintf(ref_outfolder, "%s/%s", argv[3], outputfolder);
                sprintf(opt_outfolder, "%s/%s", argv[4], outputfolder);

                sprintf(ref_shotfolder, "%s/shot.%2.1f.%05d", ref_outfolder, waveletFreq, shot);
                sprintf(opt_shotfolder, "%s/shot.%2.1f.%05d", opt_outfolder, waveletFreq, shot);

                kernel( RTM_KERNEL, waveletFreq, shot, 
                        ref_outfolder, opt_outfolder,
                        ref_shotfolder, ref_shotfolder );

                fprintf(stderr, "       %d-th shot processed\n", shot);
            }

            //TODO: in gather_shots instead of 'gathering' all shots, just compare accumulated solutions!
            //gather_shots( outputfolder, waveletFreq, nshots, numberOfCells );

            for(int test=0; test<ntest; test++)
            {
                fprintf(stderr, "Processing %d-th test iteration.\n", test);
                for(int shot=0; shot<nshots; shot++)
                {
                    char ref_shotfolder[200], opt_shotfolder[200];
                    char ref_outfolder[200],  opt_outfolder[200];
                    sprintf(ref_outfolder, "%s/%s", argv[3], outputfolder);
                    sprintf(opt_outfolder, "%s/%s", argv[4], outputfolder);

                    sprintf(ref_shotfolder, "%s/test.%05d.shot.%2.1f.%05d", ref_outfolder, test, waveletFreq, shot);
                    sprintf(opt_shotfolder, "%s/test.%05d.shot.%2.1f.%05d", opt_outfolder, test, waveletFreq, shot);

                    kernel( FM_KERNEL , waveletFreq, shot, 
                            ref_outfolder, opt_outfolder,
                            ref_shotfolder, opt_shotfolder );
                
                    fprintf(stderr, "       %d-th shot processed\n", shot);
                }
            }
        } /* end of grad loop */

    } /* end of frequency loop */


    fprintf(stderr, "-------- FWI program Finished ------------------- \n");

    return 0;
}

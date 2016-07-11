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

void kernel( propagator_t propagator, real waveletFreq, int shotid, char* outputfolder)
{
    int stacki;
    real dt,dz,dx,dy;
    integer dimmz, dimmx, dimmy;
    int forw_steps, back_steps;

    char shotfolder[200];
    sprintf( shotfolder, "%s/shot.%05d", outputfolder, shotid);
    load_shot_parameters( shotid, &stacki, &dt, &forw_steps, &back_steps, &dz, &dx, &dy, &dimmz, &dimmx, &dimmy, outputfolder);

    const integer numberOfCells = dimmz * dimmx * dimmy;

    /* set LOCAL integration limits */
    const integer nz0 = 0;
    const integer ny0 = 0;
    const integer nx0 = 0;
    const integer nzf = dimmz;
    const integer nxf = dimmx;
    const integer nyf = dimmy;

		real    *rho;
    v_t     v;
    s_t     s;
    coeff_t coeffs;

		fprintf(stderr, "Length of local arrays is " I " cells\n", numberOfCells);

		/* allocate shot memory */
    alloc_memory_shot  ( numberOfCells, &coeffs, &s, &v, &rho);

		/* load initial model from a binary file */
		load_initial_model ( waveletFreq, numberOfCells, &coeffs, &s, &v, rho);


		/* Allocate memory for IO buffer */
    real* io_buffer = (real*) __malloc( ALIGN_REAL, numberOfCells * sizeof(real) * WRITTEN_FIELDS );

    /* inspects every array positions for leaks. Enabled when DEBUG flag is defined */
    check_memory_shot  ( numberOfCells, &coeffs, &s, &v, rho);

		/* some variables for timming */
		double start_t, end_t;

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
                        numberOfCells,
                        dimmz, dimmx);

				end_t = dtime();

        fprintf(stderr, "Forward propagation finished in %lf seconds\n", \
												end_t - start_t );

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

        fprintf(stderr, "Backward propagation finished in %lf seconds\n", \
												end_t - start_t );

        #ifdef DO_NOT_PERFORM_IO
          fprintf(stderr, "Warning: we are not doing any IO here (%s)\n", __FUNCTION__ );
      	  
				#else

          char fnameGradient[300];
  				char fnamePrecond[300];
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

        fprintf(stderr, "Forward Modelling finished in %lf seconds\n",  \
												end_t - start_t );
       
			 	break;
      }
      default:
      {
        fprintf(stderr, "Invalid propagation identifier\n");
        abort();
      }
    }

    // liberamos la memoria alocatada en el shot
   	free_memory_shot  ( &coeffs, &s, &v, &rho);
    __free( io_buffer );
    
		fprintf(stderr, "Shot memory free'd\n");
};

void gather_shots( char* outputfolder, const int nshots, const int numberOfCells )
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
  	sprintf( readfilename, "%s/shot.%05d/precond_%05d.dat", outputfolder, shot, shot);

		fprintf(stderr, "Reading preconditioner file %s\n", readfilename );

		FILE* freadfile = safe_fopen( readfilename, "rb", __FILE__, __LINE__ );
		safe_fread ( readbuffer, sizeof(real), numberOfCells * 12, freadfile, __FILE__, __LINE__ );

		#pragma omp parallel for
		#pragma simd
		for( int i = 0; i < numberOfCells * WRITTEN_FIELDS; i++)
			sumbuffer[i] += readbuffer[i];
		fclose (freadfile);
	}

	char precondfilename[300];
  sprintf( precondfilename, "%s/Preconditioner", outputfolder);
	FILE* precondfile = safe_fopen( precondfilename, "wb", __FILE__, __LINE__ );
	safe_fwrite ( sumbuffer, sizeof(real), numberOfCells * 12, precondfile, __FILE__, __LINE__ );
	safe_fclose( precondfilename, precondfile, __FILE__, __LINE__ );

	end_t = dtime();

	fprintf(stderr, "Gatering process for preconditioner %s \ 
									 completed in: %lf seconds\n",          \
									 precondfilename, end_t - start_t  );



	/* ---------  GLOBAL GRADIENT ACCUMULATION --------- */
	fprintf(stderr, "Gathering local gradient fields\n");

	start_t = dtime();

	/* set buffer positions to zero */
	memset ( sumbuffer, 0, numberOfCells * sizeof(real) * WRITTEN_FIELDS );

	for( int shot=0; shot < nshots; shot++)
	{
		char readfilename[300];
  	sprintf( readfilename, "%s/shot.%05d/gradient_%05d.dat", outputfolder, shot, shot);

		fprintf(stderr, "Reading gradient file %s\n", readfilename );

		FILE* freadfile = safe_fopen( readfilename, "rb", __FILE__, __LINE__ );
		safe_fread ( readbuffer, sizeof(real), numberOfCells * WRITTEN_FIELDS, freadfile, __FILE__, __LINE__ );

		#pragma omp parallel for
		#pragma simd
		for( int i = 0; i < numberOfCells * WRITTEN_FIELDS; i++)
			sumbuffer[i] += readbuffer[i];

		fclose (freadfile);
	}

	char gradientfilename[300];
  sprintf( gradientfilename, "%s/Gradient", outputfolder);
	FILE* gradientfile = safe_fopen( gradientfilename, "wb", __FILE__, __LINE__ );
	safe_fwrite ( sumbuffer, sizeof(real), numberOfCells * WRITTEN_FIELDS, gradientfile, __FILE__, __LINE__ );
	safe_fclose( gradientfilename, gradientfile, __FILE__, __LINE__ );

	end_t = dtime();

	fprintf(stderr, "Gatering process for gradient %s       \ 
									 completed in: %lf seconds\n",          \
									 precondfilename, end_t - start_t  );

	__free(  sumbuffer);
	__free( readbuffer);
#endif

};

int main(int argc, const char *argv[])
{
  real lenz,lenx,leny,vmin,srclen,rcvlen;
  char outputfolder[200];

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

		fprintf(stderr, "stacki value is %d.                \
			 							 forward propagation steps: %d      \
										 backward propagation steps: %d\n", \
										 stacki, forw_steps, back_steps);

		for(int grad=0; grad<ngrads; grad++) /* iteracion de inversion */
		{
			fprintf(stderr, "Processing %d-th gradient iteration.\n", grad);

			for(int shot=0; shot<nshots; shot++)
			{
				char shotfolder[200];
				sprintf(shotfolder, "%s/shot.%05d", outputfolder, shot);
				create_folder( shotfolder );
				
				store_shot_parameters ( shot, &stacki, &dt, &forw_steps, &back_steps, 
																&dz, &dx, &dy, 
																&dimmz, &dimmx, &dimmy, 
																outputfolder);

				kernel( RTM_KERNEL, waveletFreq, shot, outputfolder);
				
				fprintf(stderr, "       %d-th shot processed\n", shot);
				// update_shot()
			}

			gather_shots( outputfolder, nshots, numberOfCells );

			for(int test=0; test<ntest; test++)
			{
				fprintf(stderr, "Processing %d-th test iteration.\n", test);
				for(int shot=0; shot<nshots; shot++)
				{
					char shotfolder[200];
					sprintf(shotfolder, "%s/shot.%05d", outputfolder, shot);
					create_folder( shotfolder );
					
					store_shot_parameters ( shot, &stacki, &dt, &forw_steps, &back_steps, 
																	&dz, &dx, &dy, 
																	&dimmz, &dimmx, &dimmy, 
																	outputfolder);

					kernel( FM_KERNEL , waveletFreq, shot, outputfolder);
				
					fprintf(stderr, "       %d-th shot processed\n", shot);
				}
			}
		} /* end of test loop */

  } /* end of frequency loop */


  fprintf(stderr, "-------- FWI program Finished ------------------- \n");
  return 0;
}

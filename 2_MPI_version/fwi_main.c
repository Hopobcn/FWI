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
 * Para generar el fichero de la fuente tengo que emplear las funciones del fichero
 * /home/srodrigu/Remote/scratch/BSITLocal/trunk/system/support/bscgeo/src/wavelet.c
 */

void kernel( propagator_t propagator, real waveletFreq, int shotid, char* outputfolder)
{   
	/* find ourselves into the MPI space */
	int localRank, Subdomains;
	MPI_Comm_size( MPI_COMM_WORLD, &Subdomains);
	MPI_Comm_rank( MPI_COMM_WORLD, &localRank);
  
	/* local variables */
	int stacki;
  real dt,dz,dx,dy;
  integer dimmz, dimmx, dimmy;
  int forw_steps, back_steps;
  
	/* load simulation configuration parameters */	
  char shotfolder[200];
  sprintf( shotfolder, "%s/shot.%05d", outputfolder, shotid);
  load_shot_parameters( shotid, &stacki, &dt, &forw_steps, &back_steps, &dz, &dx, &dy, &dimmz, &dimmx, &dimmy, outputfolder);


	/* number of cell for each MPI sub-problem */
	const integer planesPerSubdomain = (dimmy/Subdomains);
	const integer numberOfCells      = (dimmz + 2*HALO) * (dimmx + 2*HALO) * ( planesPerSubdomain + 2*HALO);
  
	real    *rho;
  v_t     v;
  s_t     s;
  coeff_t coeffs;

	log_info( "Computing " I " planes, and " I " cells", planesPerSubdomain, numberOfCells );
  
	/* allocate shot memory */
  alloc_memory_shot  ( numberOfCells, &coeffs, &s, &v, &rho);
  
	/* load initial model from a binary file */
	load_initial_model ( waveletFreq, dimmz, dimmx, dimmy, &coeffs, &s, &v, rho);
  
  /* reservamos memoria para el buffer de I/O */
  real* io_buffer = (real*) __malloc( ALIGN_REAL, numberOfCells * sizeof(real) * WRITTEN_FIELDS );

	/* Compute local computing limits */
	const integer y0 = planesPerSubdomain * localRank; // (dimmz + 2*HALO) * (dimmx + 2*HALO) * planesPerSubdomain *  localRank;
	const integer yF = y0 + planesPerSubdomain; // (dimmz + 2*HALO) * (dimmx + 2*HALO) * planesPerSubdomain * (localRank +1);

	log_info ( "We have to compute from y=" I " to=" I ". Planes per subdomain " I " ", y0, yF, planesPerSubdomain);


	log_info ("Starting wave propagation");
  
	switch( propagator )
  {
    case( RTM_KERNEL ):
    {
			double time_forward_b, time_forward_e, time_backward_b, time_backward_e;

			time_forward_b = dtime();
      propagate_shot ( FORWARD, 
                      v, s, coeffs, rho,
                      forw_steps, back_steps -1, // ?? limits ??
                      dt,dz,dx,dy,
                      HALO, dimmz + HALO, HALO, dimmx + HALO, y0, yF,
                      stacki,
                      shotfolder,
                      io_buffer,
                      numberOfCells,
                      dimmz, dimmx);
			
			time_forward_e = dtime();

      log_info ( "FORWARD propagation completed in %lf seconds", time_forward_e - time_forward_b);
			
			time_backward_b = dtime();
      propagate_shot ( BACKWARD, 
                      v, s, coeffs, rho,
                      forw_steps, back_steps -1, // ?? domain limits ?? 
                      dt,dz,dx,dy,
                      HALO, dimmz + HALO, HALO, dimmx + HALO, y0, yF,
                      stacki,
                      shotfolder,
                      io_buffer,
                      numberOfCells,
                      dimmz, dimmx);

			time_backward_e = dtime();

      log_info ( "BACKWARD propagation completed in %lf seconds", time_backward_e - time_backward_b);
/*
      char fnameGradient[300];
			char fnamePrecond[300];
      sprintf( fnameGradient, "%s/gradient_%05d.dat", shotfolder, shotid );
      sprintf( fnamePrecond , "%s/precond_%05d.dat" , shotfolder, shotid );

      FILE* fgradient = safe_fopen( fnameGradient, "wb", __FILE__, __LINE__ );
      FILE* fprecond  = safe_fopen( fnamePrecond , "wb", __FILE__, __LINE__ );

			fprintf(stderr, "Guardando el gradiente en %s\n", fnameGradient );
      safe_fwrite( io_buffer, sizeof(real), numberOfCells * WRITTEN_FIELDS, fgradient, __FILE__, __LINE__ );
			
			fprintf(stderr, "Guardando el precondicionador en %s\n", fnamePrecond);
      safe_fwrite( io_buffer, sizeof(real), numberOfCells * WRITTEN_FIELDS, fprecond , __FILE__, __LINE__ );
      
      safe_fclose( fnameGradient, fgradient, __FILE__, __LINE__ );
      safe_fclose( fnamePrecond , fprecond , __FILE__, __LINE__ );
*/
      break;
    }
    case( FM_KERNEL  ):
    {
			double time_fm_b, time_fm_e;
		 
			time_fm_b	= dtime();
      
			propagate_shot ( FWMODEL, 
                      v, s, coeffs, rho,
                      forw_steps, back_steps -1, // ?? domain limits ??
                      dt,dz,dx,dy,
                      HALO, HALO + dimmz, HALO, HALO + dimmx, 0, dimmy,
                      stacki,
                      shotfolder,
                      io_buffer,
                      numberOfCells,
                      dimmz, dimmx);

			time_fm_e = dtime();
      
			log_info ( "FMODELLING propagation completed in %lf seconds", time_fm_e - time_fm_b );
      break;
    }
    default:
    {
      log_error ( "Invalid propagation identifier" );
      abort();
    }
  }

    // liberamos la memoria alocatada en el shot
   	free_memory_shot  ( &coeffs, &s, &v, &rho);
    __free( io_buffer );
};

void gather_shots( char* outputfolder, const int nshots, const int numberOfCells )
{
	log_info("Gathering shots...");

#ifdef DO_NO_PERFORM_IO
	/* ---------  GLOBAL PRECONDITIONER ACCUMULATION --------- */
	log_info ("Gathering all shots to generate the preconditioner file for the next iteration");


	real* sumbuffer  = (real*)  __malloc( ALIGN_REAL, numberOfCells * sizeof(real) * WRITTEN_FIELDS ); /* TODO */
	real* readbuffer = (real*)  __malloc( ALIGN_REAL, numberOfCells * sizeof(real) * WRITTEN_FIELDS ); /* TODO */

	/* set buffer positions to zero */
	memset ( sumbuffer, 0, numberOfCells * sizeof(real) * WRITTEN_FIELDS );

	for( int shot=0; shot < nshots; shot++)
	{
		char readfilename[300];
  	sprintf( readfilename, "%s/shot.%05d/precond_%05d.dat", outputfolder, shot, shot);

		log_info ( "Reading partial preconditioner file %s", readfilename );

		FILE* freadfile = safe_fopen( readfilename, "rb", __FILE__, __LINE__ );
		safe_fread ( readbuffer, sizeof(real), numberOfCells * WRITTEN_FIELDS, freadfile, __FILE__, __LINE__ );

		#pragma omp parallel for
		#pragma simd
		for( int i = 0; i < numberOfCells * WRITTEN_FIELDS; i++)
			sumbuffer[i] += readbuffer[i];
		fclose (freadfile);
	}

	char precondfilename[300];
  sprintf( precondfilename, "%s/Preconditioner", outputfolder);
	FILE* precondfile = safe_fopen( precondfilename, "wb", __FILE__, __LINE__ );
	safe_fwrite ( sumbuffer, sizeof(real), numberOfCells * WRITTEN_FIELDS, precondfile, __FILE__, __LINE__ );
	safe_fclose( precondfilename, precondfile, __FILE__, __LINE__ );

	log_info ("Preconditioner file '%s' generated successfully", precondfilename );



	/* ---------  GLOBAL GRADIENT ACCUMULATION --------- */
	log_info ("Gathering all shots to generate the gradient file for the next iteration");

	/* set buffer positions to zero */
	memset ( sumbuffer, 0, numberOfCells * sizeof(real) * WRITTEN_FIELDS );

	for( int shot=0; shot < nshots; shot++)
	{
		char readfilename[300];
  	sprintf( readfilename, "%s/shot.%05d/gradient_%05d.dat", outputfolder, shot, shot);

		log_info ( "Reading partial gradient file %s", readfilename );

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

	log_info ("Gradient file '%s' generated successfully", gradientfilename );

	__free(  sumbuffer);
	__free( readbuffer);

#endif /* perform io conditional compilation */

	log_info("Shots gathering completed successfully");
};

int main(int argc, char* argv[] )
{
	MPI_Init ( &argc, &argv );
	int mpi_rank, subdomains;
	MPI_Comm_size( MPI_COMM_WORLD, &subdomains);
	MPI_Comm_rank( MPI_COMM_WORLD, &mpi_rank);

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
		real waveletFreq = frequencies[i];
		
		/* compute discretization deltas, 16 == puntos por longitud de onda */
		real dx = vmin / (16.0 * waveletFreq);
		real dy = vmin / (16.0 * waveletFreq);
		real dz = vmin / (16.0 * waveletFreq);
		
		/* number of cells along axis */
		integer dimmz = ceil( lenz / dz );
		integer dimmy = ceil( leny / dy );
		integer dimmx = ceil( lenx / dx );
		
		log_info ( "Domain dimensions (dimm x, y, z) " I I I ", delta of space (dx,dy,dz) %f %f %f vmin %f", dimmx, dimmy, dimmz, dx, dy, dz, vmin);
		
		/* compute delta T */
		real dt = 68e-6 * dx;
		
		/* dynamic I/O */
		int stacki = floor(  0.25 / (2.5 * waveletFreq * dt) );
		
		const integer numberOfCells = (dimmy +2*HALO) * (dimmx +2*HALO) * (dimmx +2*HALO);
		const integer VolumeMemory  = numberOfCells * sizeof(real) * 58 * WRITTEN_FIELDS;
		log_info ( "Size of the volume is " I " bytes (%f GB)", mpi_rank, VolumeMemory, TOGB(VolumeMemory) );
		
		/* compute time steps */
		int forw_steps = max_int ( IT_FACTOR * (srclen/dt), 1);
		int back_steps = max_int ( IT_FACTOR * (rcvlen/dt), 1);
		
		log_info ( "stacki = %d,  number of forward and backward iterations (%d,%d)", stacki, forw_steps, back_steps);
		
		forw_steps = 2;
		back_steps = 2;
		
		for(int grad=0; grad<ngrads; grad++) /* iteracion de inversion */
		{
			log_info ( "Processing gradient loop (iter %d)", grad);
			
			for(int shot=0; shot<nshots; shot++)
			{
				if ( mpi_rank == 0 ) 
				{
					char shotfolder[200];
					sprintf(shotfolder, "%s/shot.%05d", outputfolder, shot);
					create_folder( shotfolder );
					store_shot_parameters( shot, &stacki, &dt, &forw_steps, &back_steps, &dz, &dx, &dy, &dimmz, &dimmx, &dimmy, outputfolder);
				}

				MPI_Barrier( MPI_COMM_WORLD );

				kernel( RTM_KERNEL, waveletFreq, shot, outputfolder);

				// update_shot()
				
				log_info ( "Gradient iterations for shot id %d at %.2f Hz completed", shot, waveletFreq );
			}
		
			// apilamos los shots
			MPI_Barrier( MPI_COMM_WORLD );
			
			if ( mpi_rank == 0 ) { 
				gather_shots( outputfolder, nshots, numberOfCells );	
			}
			
			MPI_Barrier( MPI_COMM_WORLD );
			
			for(int test=0; test<ntest; test++)
			{
				for(int shot=0; shot<nshots; shot++)
				{
					if ( mpi_rank == 0)
					{
						char shotfolder[200];
						sprintf(shotfolder, "%s/shot.%05d", outputfolder, shot);
						create_folder( shotfolder );	
						store_shot_parameters( shot, &stacki, &dt, &forw_steps, &back_steps, &dz, &dx, &dy, &dimmz, &dimmx, &dimmy, outputfolder);
					}

					MPI_Barrier( MPI_COMM_WORLD );
					kernel( FM_KERNEL , waveletFreq, shot, outputfolder);
					
					log_info ( "Forward modelling iterations for shot id %d at %.2f Hz completed", shot, waveletFreq );
				}
			}
		} /* end of test loop */
		
	} /* end of frequency loop */

	MPI_Barrier(MPI_COMM_WORLD);

	log_info ( "-------------------------- Program Finished ------------------------ \n");
	MPI_Finalize();

  return 0;
}


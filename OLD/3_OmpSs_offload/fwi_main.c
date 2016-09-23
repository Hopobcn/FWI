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
#include "fwi_sched.h"

/*
 * In order to generate a source for injection,
 * /system/support/bscgeo/src/wavelet.c
 * functions can be used.
 */




/*
 * propagator: kernel mode
 * waveletFreq: wavelet frequency in Hz.
 * shotid: shot identificator (integer)
 * outputfolder:
 * nworkers: number of workers required to compute the shot.
 * ppw: number of y-planes per worker.
 */
void kernel ( propagator_t propagator, 
							real waveletFreq, 
							int shotid, 
							char* outputfolder, 
							int nworkers, 
							int ppw )
{

#ifdef DEBUG
	fprintf(stderr, "%s input arguments: %d %f %d %s %d %d\n", __FUNCTION__, propagator, waveletFreq, shotid, outputfolder, nworkers, ppw );
#endif


	/* allocate slave nodes */
		booster_alloc_t workers = allocate_workers( nworkers, 1, shotid );

 		/* load simulation parameters */		
    real dt,dz,dx,dy;
    integer dimmz, dimmx, dimmy;
    int stacki, forw_steps, back_steps;

    char shotfolder[200];
    sprintf( shotfolder, "%s/shot.%05d", outputfolder, shotid);
    load_shot_parameters( shotid, &stacki, &dt, 
													&forw_steps, &back_steps, 
													&dz, &dx, &dy, 
													&dimmz, &dimmx, &dimmy, 
													outputfolder);

		for(int worker = 0; worker < nworkers; worker++)
		{
			/* compute the number of planes that needs to be processed */
			if ( worker < (nworkers -1) || (dimmy % ppw) == 0)
				dimmy = ppw;
			else
				dimmy = dimmy % ppw;

			fprintf(stderr, "%d-th worker is processing %d y-planes\n", worker, dimmy);
    	
			const integer numberOfCells = dimmz * dimmx * dimmy;
			
			#pragma omp task onto(workers.intercomm, worker) in(propagator, shotid, [200]shotfolder) copy_deps
			{
				real    *rho;
    		v_t     v;
    		s_t     s;
    		coeff_t coeffs;

				/* allocate shot memory */
    		alloc_memory_shot  ( numberOfCells, &coeffs, &s, &v, &rho);

				/* load initial model from a binary file */
				load_initial_model ( waveletFreq, numberOfCells, &coeffs, &s, &v, rho);

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
                        0, dimmz, 0, dimmx, 0, dimmy,
                        stacki,
                        shotfolder,
                        NULL,
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
                        0, dimmz, 0, dimmx, 0, dimmy,
                        stacki,
                        shotfolder,
                        NULL,
                        numberOfCells,
                        dimmz, dimmx);

						end_t = dtime();

        		fprintf(stderr, "Backward propagation finished in %lf seconds\n", \
												end_t - start_t );

						
						/* store gradient and preconditioner fields */	
						store_field( shotfolder, shotid, GRADIENT      , &v, numberOfCells );
						store_field( shotfolder, shotid, PRECONDITIONER, &v, numberOfCells );
        		
						break;
      		}
      		case( FM_KERNEL  ):
      		{
						start_t = dtime();

        		propagate_shot ( FWMODEL,
                        v, s, coeffs, rho,
                        forw_steps, back_steps -1,
                        dt,dz,dx,dy,
                        0, dimmz, 0, dimmx, 0, dimmy,
                        stacki,
                        shotfolder,
                        NULL,
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
				/* deallocate shot memory */
   			free_memory_shot  ( &coeffs, &s, &v, &rho);
			} /* end of ompss pragma running on workers*/
		} /* end of work scheduling loop */
		#pragma omp taskwait
    
		deep_booster_free(&workers.intercomm);
};


int main(int argc, char *argv[])
{
	/* inialize nanos runtime */
	nanos_mpi_init( &argc, &argv );
	MPI_Comm spawn_comm = MPI_COMM_WORLD;

	schedule_t S = load_schedule( argv[1] );

  for(int i=0; i<S.nfreqs; i++)
  {
		fprintf(stderr, "At this frequency, we'll allocate %d slaves and %d workers\n", S.nshots, S.nworkers[i] );
		
		real freq        = S.freq[i];
		integer stacki   = S.stacki[i];
		real dt          = S.dt[i];
		integer forws    = S.forws[i];
		integer backs    = S.backs[i];
		real dz          = S.dz[i];
		real dx          = S.dx[i];
		real dy          = S.dy[i];
		integer dimmz    = S.dimmz[i];
		integer dimmx    = S.dimmx[i];
		integer dimmy    = S.dimmy[i];
		integer ppd      = S.ppd[i];
		integer nworkers = S.nworkers[i];


		/* allocate slave nodes */
		booster_alloc_t slaves = allocate_slaves( S.nshots );

		for(int grad=0; grad<S.ngrads; grad++) /* inversion iterations */
		{
			fprintf(stderr, "Processing %d-th gradient iteration.\n", grad);

			for(int shot=0; shot<S.nshots; shot++)
			{
				#pragma omp task onto(slaves.intercomm, shot) in([200]S.outputfolder) label(rtm_kernel) copy_deps
				{

					fprintf(stderr, "stacki %d forws %d backs %d\n", stacki, forws, backs );

					char shotfolder[200];
					sprintf( shotfolder, "%s/shot.%05d", S.outputfolder, shot);
					create_folder( shotfolder );
					
					store_shot_parameters ( shot, &stacki, &dt, &forws, &backs, 
																	&dz, &dx, &dy, 
																	&dimmz, &dimmx, &dimmy, 
																	S.outputfolder);

					kernel( RTM_KERNEL, freq, shot, S.outputfolder, nworkers, ppd );
					
					fprintf(stderr, "       %d-th shot processed\n", shot);
				}
			}
			#pragma omp taskwait

			/* shot gathering */
			// gather_shots( outputfolder, nshots, numberOfCells );
	
			for(int test=0; test<S.ntests; test++)
			{
				fprintf(stderr, "Processing %d-th test iteration.\n", S.ntests);
				 
				for(int shot=0; shot<S.nshots; shot++)
				{
					#pragma omp task onto(slaves.intercomm, shot) in([200]S.outputfolder) label(test_iterations) copy_deps
					{
						char shotfolder[200];
						sprintf(shotfolder, "%s/shot.%05d", S.outputfolder, shot);
						create_folder( shotfolder );
					
						store_shot_parameters ( shot, &stacki, &dt, &forws, &backs, 
																	&dz, &dx, &dy, 
																	&dimmz, &dimmx, &dimmy, 
																	S.outputfolder);

						kernel( FM_KERNEL , freq, shot, S.outputfolder, nworkers, ppd );
				
						fprintf(stderr, "       %d-th shot processed\n", shot);
					}// end of ompss pragma
				}
				#pragma omp taskwait
			}
		} /* end of test loop */
    deep_booster_free(&slaves.intercomm);
  } /* end of frequency loop */
  fprintf(stderr, "-------- FWI propagator Finished ------------------- \n");

	/* clean UP */
	int rank, spawn_rank, size;
	MPI_Comm_rank( MPI_COMM_WORLD, &rank );
	MPI_Comm_size( MPI_COMM_WORLD, &size );
	MPI_Comm_rank( spawn_comm, &spawn_rank );
	
	/* wait for all the processes to finish */
	fprintf(stderr, "Waiting for all processes to finish\n");
	MPI_Barrier( MPI_COMM_WORLD );

	fprintf(stderr, "Waiting for nanos runtime to finish\n");
	/* explicitely finalize nanos runtime */
	nanos_mpi_finalize();

	fprintf(stderr, "End of the program.\n");
	
  return 0;
}

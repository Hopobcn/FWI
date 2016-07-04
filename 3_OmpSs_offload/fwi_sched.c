/*
 * =====================================================================================
 *
 *       Filename:  fwi_sched.c
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  04/04/2016 09:03:09 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Samuel Rodriguez Bernabeu (), samuel.rodriguez@bsc.es
 *   Organization:  Barcelona Supercomputing Center (BSC)
 *
 * =====================================================================================
 */


#include "fwi_sched.h"

void schedule_free( schedule_t S )
{
	free ( S.freq );
	free ( S.forws );
	free ( S.backs );
	free ( S.stacki );
	free ( S.dt );
	free ( S.dz );
	free ( S.dx );
	free ( S.dy );
	free ( S.dimmz );
	free ( S.dimmx );
	free ( S.dimmy );
	free ( S.ppd );
	free ( S.nworkers );
};

schedule_t load_schedule( const char* filename ) 
{
	/* open file to store the schedule */
	char* fwipath = read_env_variable("FWIDIR");
	int pathlen = strlen(fwipath) + 200;
	char filepath[pathlen];
	sprintf(filepath, "%s/SetupParams/%s", fwipath, filename);

	FILE* fschedule = fopen( filepath, "r");

	if ( fschedule == NULL ){
		fprintf(stderr, "Cant open scheduling file %s\n", filepath );
		abort();
	}

	schedule_t S;

	/* read the number of frequencies to be processed */
	fscanf( fschedule, "%d\n", &S.nfreqs );
	fprintf(stderr, "Number of frequencies %d\n", S.nfreqs);

	/* read the number of shots to be processed */
	fscanf( fschedule, "%d\n", &S.nshots );
	fprintf(stderr, "Number of shots %d\n", S.nshots);
	
	/* read the number of gradient iterations */
	fscanf( fschedule, "%d\n", &S.ngrads );
	fprintf(stderr, "Number of gradient iterations %d\n", S.ngrads);

	/* read the number of test iterations */
	fscanf( fschedule, "%d\n", &S.ntests );
	fprintf(stderr, "Number of test iterations %d\n", S.ntests);

	/* read the name of the output directory */
	fscanf( fschedule, "%s\n", &S.outputfolder );

	fprintf(stderr, "Output directory path: %s\n", S.outputfolder );

	/* allocate memory for the rest of the parameters */
	S.freq     = (real*   ) malloc( S.nfreqs * sizeof(real   ));
	S.forws    = (integer*) malloc( S.nfreqs * sizeof(integer));
	S.backs    = (integer*) malloc( S.nfreqs * sizeof(integer));
	S.stacki   = (integer*) malloc( S.nfreqs * sizeof(integer));
	S.dt       = (real*   ) malloc( S.nfreqs * sizeof(real   ));
	S.dz       = (real*   ) malloc( S.nfreqs * sizeof(real   ));
	S.dx       = (real*   ) malloc( S.nfreqs * sizeof(real   ));
	S.dy       = (real*   ) malloc( S.nfreqs * sizeof(real   ));
	S.dimmz    = (integer*) malloc( S.nfreqs * sizeof(integer));
	S.dimmx    = (integer*) malloc( S.nfreqs * sizeof(integer));
	S.dimmy    = (integer*) malloc( S.nfreqs * sizeof(integer));
	S.ppd      = (integer*) malloc( S.nfreqs * sizeof(integer));
	S.nworkers = (integer*) malloc( S.nfreqs * sizeof(integer));

	/* read the list of simulation parameters */
	for( int i=0; i < S.nfreqs; i++ ) 
	{
		fscanf( fschedule, "%f %d %d %d %f %f %f %f %d %d %d %d %d\n", 
				&S.freq[i], &S.forws[i], &S.backs[i], &S.stacki[i], 
				&S.dt[i], &S.dz[i], &S.dy[i], &S.dx[i], 
				&S.dimmz[i], &S.dimmx[i], &S.dimmy[i], 
				&S.ppd[i], &S.nworkers[i]);
	}

	/* show up a summary */
	fprintf(stderr, "Info: number of frequencies %d, number of shots %d\n", S.nfreqs, S.nshots);

	for( int i=0; i < S.nfreqs; i++ ) 
	{
		fprintf( stderr, "%f %d %d %d %f %f %f %f %d %d %d %d %d\n", 
				S.freq[i], S.forws[i], S.backs[i], S.stacki[i], 
				S.dt[i], S.dz[i], S.dy[i], S.dx[i], 
				S.dimmz[i], S.dimmx[i], S.dimmy[i], 
				S.ppd[i], S.nworkers[i]);
	}
	/* clean up and resume */
	fclose( fschedule ); fschedule = NULL;

	return (S);
};


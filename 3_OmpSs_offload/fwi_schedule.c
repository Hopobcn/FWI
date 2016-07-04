/*
 * =====================================================================================
 *
 *       Filename:  fwi_estimator.c
 *
 *    Description:  Estimates the number of workers (accelators) needed to carry
 *    							on the simulation. To that end, we use the max. wavelet
 *    							frequency of the simulation and the available memory on the
 *    							accelerator device (in GB).
 *
 *        Version:  1.0
 *        Created:  10/12/15 10:33:40
 *       Revision:  none
 *       Compiler:  icc/gcc
 *
 *         Author:  Samuel Rodriguez Bernabeu, samuel.rodriguez(at)bsc.es
 *   Organization:  Barcelona Supercomputing Center
 *
 * =====================================================================================
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <math.h>

typedef float real;
typedef int   integer;

const integer HALO = 4;
const real IT_FACTOR = 0.02;

char* read_env_variable (const char* varname);

void read_fwi_parameters (char *fname,
                          real *lenz,
                          real *lenx,
                          real *leny,
                          real *vmin,
                          real *srclen,
                          real *rcvlen,
													int  *nshots,
													int  *ngrads,
													int  *ntests,
													int  *workermem,
													int  *slavemem,
                          char *outputfolder);

void load_freqlist( char* filename, int *nfreqs, real **freqlist );

integer roundup(integer number, integer multiple);

int max_int( int a, int b);

typedef struct{
	integer nfreqs;
	integer nshots;
	integer ngrads;
	integer ntests;
	char    outputfolder[200];

	real    *freq;
	integer *forws;
	integer *backs;
	integer *stacki;
	real    *dt;
	real    *dz;
	real    *dx;
	real    *dy;	
	integer    *dimmz;
	integer    *dimmx;
	integer    *dimmy;
	integer *ppd;
	integer *nworkers;
} schedule_t;

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
}

schedule_t load_schedule( void ) 
{
	/* open file to store the schedule */
	FILE* fschedule = fopen("../SetupParams/fwi_schedule.txt", "r");

	if ( fschedule == NULL ){
		fprintf(stderr, "Cant create scheduling file\n");
		abort();
	}

	schedule_t S;

	/* read the number of frequencies to be processed */
	fscanf( fschedule, "%d\n", &S.nfreqs );
	fprintf(stderr, "Number of different frequencies %d\n", S.nfreqs);

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


	fprintf(stderr, "We have to read %d rows of the file\n", S.nfreqs );

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





int main(int argc, char *argv[])
{
	/* Load parameters that define the simulation*/
	int nshots, ngrads, ntests, nfreqs, slavemem, workermem, nworkers;
  real lenz,lenx,leny,vmin,srclen,rcvlen,*frequencies;
	size_t device_mem;
	char outputfolder[1024];

	/* modify input parameters according to FWIDIR value */
	char* fwipath = read_env_variable("FWIDIR");
	int pathlen = strlen(fwipath) + 100;
	char paramsfile[pathlen];
	char freqsfile[pathlen];
	char schedfile[pathlen];

	sprintf( paramsfile, "%s/SetupParams/%s", fwipath, argv[1]);
	sprintf( freqsfile,  "%s/SetupParams/%s", fwipath, argv[2]);
	sprintf( schedfile,  "%s/SetupParams/fwi_schedule.txt", fwipath );


  read_fwi_parameters ( paramsfile, 
												&lenz, &lenx, &leny, &vmin, 
												&srclen, &rcvlen, 
												&nshots, &ngrads, &ntests,
												&slavemem, &workermem,
												outputfolder);

	load_freqlist( freqsfile, &nfreqs, &frequencies );

	/* open file to store the schedule */

	FILE* fschedule = fopen( schedfile, "w");

	if ( fschedule == NULL ){
		fprintf(stderr, "Cant create scheduling file\n");
		abort();
	}

	/* store the number of frequencies and shots on the schedule file */
	fprintf( fschedule, "%d\n", nfreqs );
	fprintf( fschedule, "%d\n", nshots );
	fprintf( fschedule, "%d\n", ngrads );
	fprintf( fschedule, "%d\n", ntests );
	fprintf( fschedule, "%s\n", outputfolder );

	for( int freq = 0; freq < nfreqs; freq++)
	{
		/* our calculations are based on the max frequency */
		real waveletFreq = frequencies[ freq ];
		fprintf(stderr, "Estimating resources for %f Hz...\n", waveletFreq);

		/* Deltas of space, 16 grid point per Hz */
		real dx = vmin / (16.0 * waveletFreq);
		real dy = vmin / (16.0 * waveletFreq);
		real dz = vmin / (16.0 * waveletFreq);

		/* number of cells along axis, adding HALO planes */
		integer dimmz = roundup(ceil( lenz / dz ) + 2*HALO, HALO);
		integer dimmy = roundup(ceil( leny / dy ) + 2*HALO, HALO);
		integer dimmx = roundup(ceil( lenx / dx ) + 2*HALO, HALO);

		/* compute delta of t */
		real dt = 68e-6 * dx;

		/* dynamic IO parameter */
		int stacki = floor(  0.25 / (2.5 * waveletFreq * dt) );
		
		/* compute time steps */
		int forw_steps = max_int ( IT_FACTOR * (srclen/dt), 1);
		int back_steps = max_int ( IT_FACTOR * (rcvlen/dt), 1);

#ifdef DEBUG
		fprintf(stderr, "Value of stacki %d forward steps %d backward steps %d dt %f\n", stacki, forw_steps, back_steps, dt );
		fprintf(stderr, "Value of dz = %f dx = %f dy = %f\n", dz, dx, dy);
#endif

		integer yplane_mem = dimmx * dimmz * 58 * sizeof(float); 


		/* if the accelerator memory is 0 means that it wont be used */
		if ( workermem )
		{
			device_mem = (0.8 * workermem) * 1024 * 1024 * 1024;
		}
		else
		{
			device_mem = (0.8 * slavemem ) * 1024 * 1024 * 1024; 
		}

		/* compute the number of y-planes fitting into the device */
		integer ppd = device_mem / yplane_mem;

		/* minum number of planes have to be satisfied */
		if ( ppd < 4 * HALO ){
			fprintf(stderr, "ERROR: At least 4*HALO planes must fit into a single node.\n");
			fprintf(stderr, "ERROR: You wont be able to run this configuration..\n");
		}

		/* compute how many nodes are required */
		if ( dimmy < ppd ){
			nworkers = 1;
		} else {
			nworkers = dimmy / ppd;

			if ( dimmy % ppd != 0)
				nworkers += 1;
		}

		fprintf(stderr, "\n------------------------------------------------------------------------------------\n");
		fprintf(stderr, "	There are %d y-planes to compute, each worker can hold %d of these planes\n", dimmy, ppd);
		fprintf(stderr, "	At this frequency (%f Hz) we'll need %d workers per shot\n", waveletFreq, nworkers );
		fprintf(stderr, "	There are %d shots to be computed, so %d slave nodes are needed.\n", nshots, nshots );
		fprintf(stderr, "\n------------------------------------------------------------------------------------\n");
	
		fprintf( fschedule, "%f %d %d %d %f %f %f %f %d %d %d %d %d\n", waveletFreq, forw_steps, back_steps, stacki, dt, dz, dy, dx, dimmz, dimmx, dimmy, ppd, nworkers);
	}	
	
	fclose( fschedule ); fschedule = NULL;

#ifdef DEBUG
	fprintf(stderr, "\nSchedule file will look like this...\n");
	load_schedule();
#endif

	fprintf(stderr, "\nEnd of the program.\n");
  return 0;
}

char* read_env_variable (const char* varname)
{	
	char* s = getenv(varname);
	
	if ( s == NULL )
	{
		fprintf(stderr, "%s: ERROR: unable to read  %s env. var\n", __FUNCTION__, varname);
		abort();
	}

#ifdef DEBUG
	printf("%s: %s variable value is :%s\n", __FUNCTION__, varname, s);
#endif

	return (s);
};

void read_fwi_parameters (char *fname,
                          real *lenz,
                          real *lenx,
                          real *leny,
                          real *vmin,
                          real *srclen,
                          real *rcvlen,
													int  *nshots,
													int  *ngrads,
													int  *ntests,
													int  *workermem,
													int  *slavemem,
                          char *outputfolder)
{
		fprintf(stderr, "Reading fwi parameter file %s\n", fname);

    FILE *fp = fopen(fname, "r" );

		if ( fp == NULL )
		{
			fprintf(stderr, "Cant open input parameter file %s\n", fname);
			abort();
		}


		char* buffer[200];

    fscanf( fp, "%f\n", (real*) lenz      );
    fscanf( fp, "%f\n", (real*) lenx      );
    fscanf( fp, "%f\n", (real*) leny      );
    fscanf( fp, "%f\n", (real*) vmin      );
    fscanf( fp, "%f\n", (real*) srclen    );
    fscanf( fp, "%f\n", (real*) rcvlen    );
    fscanf( fp, "%d\n", (int*)  nshots    );
    fscanf( fp, "%d\n", (int*)  ngrads    );
    fscanf( fp, "%d\n", (int*)  ntests    );
    fscanf( fp, "%d\n", (int*)  workermem );
    fscanf( fp, "%d\n", (int*)  slavemem  );
    fscanf( fp, "%s\n",  buffer           );

		/* modify output folder according to compilation flags and FWIDIR variable */

#ifdef USE_NMVE
		fprintf(stderr, "Using MNE!!\n");
		sprintf(outputfolder, "%s", "/nvme/tmp");
#else
		char* fwipath = read_env_variable("FWIDIR");
		sprintf(outputfolder, "%s/%s", fwipath, buffer);
#endif

		fprintf(stderr, "Local IO path is %s\n", outputfolder);

    fclose(fp);
};

void load_freqlist( char* filename, int *nfreqs, real **freqlist )
{
	int count  = 0;
	real freq;

	fprintf(stderr, "Reading frequency file %s\n", filename);

	FILE *freqfile = fopen( filename, "r");

	if ( freqfile == NULL )
	{
		fprintf(stderr, "Cant open frequency list file %s\n", filename );
		abort();
	}


	while( 1 )
	{
		int n = fscanf( freqfile, "%f", &freq);

		if ( n == 1 )
		{
			count += 1;
		}
	 	else if (errno != 0)
		{
			fprintf(stderr, "Error while reading freqlist file\n");
			break;
		}
		else if ( n == EOF )
		{
			break;
		}
	}


	/* Allocate memory for frequencies */
	*freqlist = (real*) malloc( count * sizeof(real));

	/* return to initial position */
	fseek( freqfile, 0, SEEK_SET);
	count = 0;



	/* read again the file, storing the wavelet frequencies this time */
	while( 1 )
	{
		int n = fscanf( freqfile, "%f", &freq);

		if ( n == 1 )
		{
			(*freqlist)[count++] = freq;
		}
	 	else if (errno != 0)
		{
			fprintf(stderr, "Error while reading freqlist file\n");
			break;
		}
		else if ( n == EOF )
		{
			break;
		}
	}
	fclose( freqfile );

	*nfreqs = count;

	fprintf(stderr, "\nNumber of detected frequencies, %d)\n", *nfreqs );
	for( int i=0; i<count; i++)
		fprintf(stderr, "     %.2f Hz\n", (*freqlist)[i] );

	fprintf(stderr, "\n");
};

integer roundup(integer number, integer multiple)
{
    if (multiple == 0)
        return number;

    int remainder = number % multiple;
    if (remainder == 0)
        return number;

    return number + multiple - remainder;
};

int max_int( int a, int b)
{
    return ((a >= b) ? a : b);
};


/*
 * =====================================================================================
 *
 *       Filename:  fwi_common.c
 *
 *    Description:
 *
 *        Version:  1.0
 *        Created:  10/12/15 10:38:45
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  YOUR NAME (),
 *   Organization:
 *
 * =====================================================================================
 */

#include "fwi_common.h"

/* extern variables declared in the header file */
const integer  WRITTEN_FIELDS =   12; /* >= 12.  */
const integer  HALO           =    4; /* >= 4    */ 
const integer  SIMD_LENGTH    =    8; /* # of real elements fitting into regs */
const real     IT_FACTOR      = 0.02;
const real     IO_CHUNK_SIZE  = 1024.f * 1024.f;

const size_t ALIGN_INT     = 16;
const size_t ALIGN_INTEGER = 16;
const size_t ALIGN_REAL    = 64;

int max_int( int a, int b)
{
    return ((a >= b) ? a : b);
};

double dtime(void)
{
	double tseconds = 0.0;
	struct timeval mytime;
	gettimeofday( &mytime, (struct timezone*) 0);
	tseconds = (double) (mytime.tv_sec + mytime.tv_usec * 1.0e-6);
	return (tseconds);
};

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

void read_fwi_parameters (const char *fname,
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
    FILE *fp = safe_fopen(fname, "r", __FILE__, __LINE__ );

    fscanf( fp, "%f\n", (real*) lenz   );
    fscanf( fp, "%f\n", (real*) lenx   );
    fscanf( fp, "%f\n", (real*) leny   );
    fscanf( fp, "%f\n", (real*) vmin   );
    fscanf( fp, "%f\n", (real*) srclen );
    fscanf( fp, "%f\n", (real*) rcvlen );
    fscanf( fp, "%d\n", (int*)  nshots );
    fscanf( fp, "%d\n", (int*)  ngrads );
    fscanf( fp, "%d\n", (int*)  ntests );
    fscanf( fp, "%d\n", (int*)  workermem );
    fscanf( fp, "%d\n", (int*)  slavemem );
    fscanf( fp, "%s\n",  outputfolder  );

    fclose(fp);
};

/*
  This function is intended to round up a number (number) to the nearest multiple of the register
  size. In this way, we assure that the dimensions of the domain are suited to the most aggressive
  compiler optimizations.
 */
integer roundup(integer number, integer multiple)
{
    if (multiple == 0)
        return number;

    int remainder = number % multiple;
    if (remainder == 0)
        return number;

    return number + multiple - remainder;
};

/*
 NAME:allocate_shot_memory
 PURPOSE: Create files to store final preconditioner and gradient results. Must be initialized with zeroes.

 outputfolder     (in) folder where snapshot data is store
 VolumeMemory     (in) memory needed to store the domain

 RETURN none
 */
void create_output_volumes(char *outputfolder, integer VolumeMemory)
{
    fprintf(stderr, "Creating output files in %s\n", outputfolder);

#ifdef DO_NOT_PERFORM_IO
    fprintf(stderr, "Warning: we are not doing any IO here (%s).\n", __FUNCTION__);
#else

		/* build full path from the environmental variable */
		char* fwipath = read_env_variable("FWIDIR");
		
		int pathlen = strlen(fwipath) + 200;

    char fnamePrecond[pathlen];
	 	char fnameGradient[pathlen];

    sprintf( fnameGradient, "%s/%s/resultGradient.res", fwipath, outputfolder);
    sprintf( fnamePrecond , "%s/%s/resultPrecond.res", fwipath, outputfolder);

		fprintf(stderr, "Creating gradient file %s\n", fnameGradient);
		fprintf(stderr, "Creating precond file %s\n",  fnamePrecond);

    FILE *fGradient = safe_fopen( fnameGradient, "wb", __FILE__, __LINE__ );
    FILE *fPrecond  = safe_fopen( fnamePrecond , "wb", __FILE__, __LINE__ );

    int numIts = ceil( VolumeMemory / IO_CHUNK_SIZE );
    fprintf(stderr, "Necesitamos realizar %d iteraciones para generar los volumenes de salida\n", numIts);

    /* create buffer array */
    real *tmparray = (real*) __malloc( ALIGN_REAL, IO_CHUNK_SIZE );

    /* perform the accumulation of the chunks */
    for (int i=0; i<numIts; i++) {
        safe_fwrite(tmparray, 1, IO_CHUNK_SIZE, fGradient, __FILE__, __LINE__ );
        safe_fwrite(tmparray, 1, IO_CHUNK_SIZE, fPrecond , __FILE__, __LINE__ );
    }

    __free(tmparray);

    // close files
    safe_fclose( fnameGradient, fGradient, __FILE__, __LINE__ );
    safe_fclose( fnamePrecond , fPrecond , __FILE__, __LINE__ );
#endif
}

/*
 NAME:create_folder
 PURPOSE:During execution creates temporal folders to organize necessary data for the execution

 folder      (in) name of the temporal folder created
 parent_rank (in) name of the rank related to the data archived in to the folder
 shotID      (in) identifier of the shot related to the data to be archived in to the folder

 RETURN none
 */
void create_folder(const char *folder)
{
    if (mkdir_p(folder) != 0) {
        fprintf(stderr,"Error  creating folder %s (%s)\n", folder, strerror(errno));
        exit(-1);
    }
    fprintf(stderr, "Folder created %s\n",folder);
};

/*
 NAME: mkdir_p
 PURPOSE: creates the hierarchy of folders requested, if they do not exist.

 RETURN 0 if successful, !=0 otherwise
 */
int mkdir_p(const char *dir)
{
    char tmp[256];
    char *p = NULL;
    size_t len;

    snprintf(tmp, sizeof(tmp),"%s",dir);
    len = strlen(tmp);

    if(tmp[len - 1] == '/')
        tmp[len - 1] = 0;

    for(p = tmp + 1; *p; p++)
        if(*p == '/') {
            *p = 0;
            int rc = mkdir(tmp, S_IRWXU);
            if (rc != 0 && errno != EEXIST) {
                fprintf(stderr,"Error creating folder %s (%s)\n", tmp, strerror(errno));
                return -1;
            }

            *p = '/';
        }

    int rc = mkdir(tmp, S_IRWXU);
    if (rc != 0 && errno != EEXIST) {
        fprintf(stderr,"Error creating folder %s (%s)\n", tmp, strerror(errno));
        return -1;
    }

    return 0;
}

void store_shot_parameters( int     shotid,
                           int     *stacki,
                           real    *dt,
                           int    *nt_fwd,
                           int    *nt_bwd,
                           real    *dz,
                           real    *dx,
                           real    *dy,
                           integer *dimmz,
                           integer *dimmx,
                           integer *dimmy,
                           char    *outputfolder)
{
		int pathlen = strlen(outputfolder) + 100;	
		char name[pathlen];

    sprintf(name, "%s/shotparams_%05d.dat",outputfolder, shotid);

    fprintf(stderr, "Storing parameters for shot %d into %s\n", shotid, name);
    FILE *fp = safe_fopen(name, "w", __FILE__, __LINE__);

    fprintf(fp, "%f\n",  (real   ) *dz     );
    fprintf(fp, "%f\n",  (real   ) *dx     );
    fprintf(fp, "%f\n",  (real   ) *dy     );
    fprintf(fp,  I"\n", (integer) *dimmz  );
    fprintf(fp,  I"\n", (integer) *dimmx  );
    fprintf(fp,  I"\n", (integer) *dimmy  );
    fprintf(fp, "%d\n",  (int    ) *nt_fwd );
    fprintf(fp, "%d\n",  (int    ) *nt_bwd );
    fprintf(fp, "%f\n",  (real   ) *dt     );
    fprintf(fp, "%d\n",  (int    ) *stacki );

    fclose(fp);
};

void load_shot_parameters( int    shotid,
                          int     *stacki,
                          real    *dt,
                          int     *nt_fwd,
                          int     *nt_bwd,
                          real    *dz,
                          real    *dx,
                          real    *dy,
                          integer *dimmz,
                          integer *dimmx,
                          integer *dimmy,
                          char    *outputfolder)
{
		int pathlen = strlen(outputfolder) + 100;
		char name[pathlen];

    sprintf(name, "%s/shotparams_%05d.dat",outputfolder, shotid);
    fprintf(stderr, "Storing parameters for shot %d into %s\n", shotid, name);

    FILE *fp = safe_fopen(name, "r", __FILE__, __LINE__);

    fscanf(fp, "%f\n",  (real*   ) dz     );
    fscanf(fp, "%f\n",  (real*   ) dx     );
    fscanf(fp, "%f\n",  (real*   ) dy     );
    fscanf(fp,  I"\n", (integer*) dimmz  );
    fscanf(fp,  I"\n", (integer*) dimmx  );
    fscanf(fp,  I"\n", (integer*) dimmy  );
    fscanf(fp, "%d\n",  (int*    ) nt_fwd );
    fscanf(fp, "%d\n",  (int*    ) nt_bwd );
    fscanf(fp, "%f\n",  (real*   ) dt     );
    fscanf(fp, "%d\n",  (int*    ) stacki );

    safe_fclose( name, fp, __FILE__, __LINE__);
};

void load_freqlist( const char* filename, int *nfreqs, real **freqlist )
{
	int count  = 0;
	real freq;

	FILE *freqfile = safe_fopen( filename, "r", __FILE__, __LINE__);

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
	*freqlist = (real*) __malloc( ALIGN_REAL, count * sizeof(real));

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

	fprintf(stderr, "\nLoaded frequencies (%d in total)\n", *nfreqs );
	for( int i=0; i<count; i++)
		fprintf(stderr, "     %.2f Hz\n", (*freqlist)[i] );

	fprintf(stderr, "\n");
};

void* __malloc( size_t alignment, const integer size)
{
  void *buffer;
  int error;

  if( (error=posix_memalign( &buffer, alignment, size)) != 0)
  {
        fprintf(stderr, "Cant allocate buffer correctly\n");
    abort();
  }

  return (buffer);
};

void __free ( void* ptr)
{
    free( ptr );
};

FILE* safe_fopen(const char *filename, char *mode, char* srcfilename, int linenumber)
{
    FILE* temp = fopen( filename, mode);

    if( temp == NULL){
        fprintf( stderr, "%s:%d Cant open filename %s, openmode '%s'\n", srcfilename, linenumber, filename, mode);
        exit(-1);
    }
    return temp;
};

void safe_fclose ( const char *filename, FILE* stream, char* srcfilename, int linenumber)
{
  if ( fclose( stream ) != 0)
  {
    fprintf(stderr, "%s:%d: Cant close file %s correctly!\n", srcfilename, linenumber, filename );
    abort();
  }

/*if ( unlink(filename)  != 0)
  {
    fprintf(stderr, "%s:%d: Cant unlink file %s correctly!\n", srcfilename, linenumber, filename );
    abort();
  }*/
};


void safe_fwrite (void *ptr, size_t size, size_t nmemb, FILE *stream, char* srcfilename, int linenumber)
{
#ifdef DO_NOT_PERFORM_IO
  fprintf(stderr, "Warning: we are not doing any IO here (%s).\n", __FUNCTION__);
#else
	size_t res = fwrite( ptr, size, nmemb, stream);

	if( res != nmemb )
	{
		fprintf(stderr, "%s:%d: Error while fwrite\n", srcfilename, linenumber );
		abort();
	}
#endif
};

void safe_fread (void *ptr, size_t size, size_t nmemb, FILE *stream, char* srcfilename, int linenumber)
{
#ifdef DO_NOT_PERFORM_IO
    fprintf(stderr, "Warning: we are not doing any IO here (%s).\n", __FUNCTION__);
#else
  if( stream == NULL ){
    fprintf(stderr, "stream is not longer valid\n");
    abort();
  }

	size_t res = fread( ptr, size, nmemb, stream);

	if( res != nmemb )
	{
		fprintf(stderr, "%s:%d: Error while fread\n", srcfilename, linenumber);
    fprintf(stderr, "Trying to read %lu elements, only %lu were recovered\n", nmemb, res);
		abort();
	}
#endif
};

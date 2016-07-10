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

/* extern variables declared in the header */
const integer  WRITTEN_FIELDS =  13;
const integer  HALO           =  8; // have to be greater than 4!
const real     IT_FACTOR      = 0.02;
const real     IO_CHUNK_SIZE  = 1024.f * 1024.f;

const size_t ALIGN_INT     = 16;
const size_t ALIGN_INTEGER = 16;
const size_t ALIGN_REAL    = 64;

extern FILE* logfile = NULL;

void log_info (const char *fmt, ...) 
{
#ifdef DEBUG
    /* locate myself into the MPI world */
    int mpi_rank;
    MPI_Comm_rank( MPI_COMM_WORLD, &mpi_rank);
    
#ifndef DEBUG_TO_STDERR
    /* build log file name */
    char logname[50];
    sprintf( logname, "mpi_%02d.log", mpi_rank);
    FILE* flog = safe_fopen( logname, "a+", __FILE__, __LINE__ );
#else
    FILE* flog = stderr;
#endif

    /* create the string from variadic input arguments */
    char str[1000];
    va_list args;
    va_start(args, fmt);
    vsprintf(str, fmt, args);
    va_end(args);
    
    /* create time string */
    char timestr[20];
    struct tm *sTm;
    time_t now = time (0);
    sTm = gmtime (&now);
    strftime ( timestr, sizeof(timestr), "%Y-%m-%d %H:%M:%S", sTm);
    
    /* print actual line to log file  */
    fprintf( flog, "%s: %s\n", timestr, str);
   
#ifndef DEBUG_TO_STDERR
    /*  close file  */
    safe_fclose( logname, flog, __FILE__, __LINE__ );
#endif
#endif
};

void log_error (const char *fmt, ...) 
{ 
    /* locate myself into the MPI world */
    int mpi_rank;
    MPI_Comm_rank( MPI_COMM_WORLD, &mpi_rank);
    
    /* build log file name */
    char logname[50];
    sprintf( logname, "mpi_%02d.log", mpi_rank);
    FILE* flog = safe_fopen( logname, "a+", __FILE__, __LINE__ );
    
    /* create the string from variadic input arguments */
    char str[1000];
    va_list args;
    va_start(args, fmt);
    vsprintf(str, fmt, args);
    va_end(args);

    /* create time string */
    char timestr[20];
    struct tm *sTm;
    time_t now = time (0);
    sTm = gmtime (&now);
    strftime ( timestr, sizeof(timestr), "%Y-%m-%d %H:%M:%S", sTm);

    /* print actual line to log file  */
    fprintf( flog, "-----> ERROR :: %s: %s\n", timestr, str);

    /*  close file  */
    safe_fclose( logname, flog, __FILE__, __LINE__ );
};




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

void read_fwi_parameters (const char *fname,
                          real *lenz,
                          real *lenx,
                          real *leny,
                          real *vmin,
                          real *srclen,
                          real *rcvlen,
                          char *outputfolder)
{   
    log_info("Loading simulation parameters");

    FILE *fp = safe_fopen(fname, "r", __FILE__, __LINE__ );
    
    CHECK( fscanf( fp, "%f\n", (real*) lenz   ) );
    CHECK( fscanf( fp, "%f\n", (real*) lenx   ) );
    CHECK( fscanf( fp, "%f\n", (real*) leny   ) );
    CHECK( fscanf( fp, "%f\n", (real*) vmin   ) );
    CHECK( fscanf( fp, "%f\n", (real*) srclen ) );
    CHECK( fscanf( fp, "%f\n", (real*) rcvlen ) );
    CHECK( fscanf( fp, "%s\n",  outputfolder  ) );
    
    fclose(fp);

    log_info("Simulation parameters loaded successfully");
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
    log_info ( "Creating output files in %s", outputfolder);

#ifndef DO_NOT_PERFOM_IO    
    char fnamePrecond[300], fnameGradient[300];
    
    sprintf( fnameGradient, "%s/resultGradient.res", outputfolder);
    sprintf( fnamePrecond , "%s/resultPrecond.res", outputfolder);
    
    FILE *fGradient = safe_fopen( fnameGradient, "wb", __FILE__, __LINE__ );
    FILE *fPrecond  = safe_fopen( fnamePrecond , "wb", __FILE__, __LINE__ );
     
    int numIts = ceil( VolumeMemory / IO_CHUNK_SIZE );
    
    /* create buffer array */
    real *tmparray = (real*) __malloc( ALIGN_INT, IO_CHUNK_SIZE );
    
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

    log_info ("Output volumes created correctly");
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
        log_error ( "Error  creating folder %s (%s)", folder, strerror(errno));
        exit(-1);
    }
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
    
    for(p = tmp + 1; *p; p++) {
        if(*p == '/') {
            *p = 0;
            int rc = mkdir(tmp, S_IRWXU);
            if (rc != 0 && errno != EEXIST) {
                fprintf(stderr,"Error creating folder %s (%s)", tmp, strerror(errno));
                return -1;
            }
            
            *p = '/';
        }
    }
    
    int rc = mkdir(tmp, S_IRWXU);
    if (rc != 0 && errno != EEXIST) {
        fprintf(stderr,"Error creating folder %s (%s)", tmp, strerror(errno));
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
    char name[200];
    
    sprintf(name, "%s/shotparams_%05d.dat",outputfolder, shotid);
    
    log_info ( "Storing parameters for shot %d into %s", shotid, name);
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
    
    safe_fclose( name,  fp, __FILE__, __LINE__ );
    
    log_info ( "Shot parameters stored correctly" );
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
    char name[200];
    
    sprintf(name, "%s/shotparams_%05d.dat",outputfolder, shotid);
    log_info ( "Storing parameters for shot %d into %s", shotid, name);

    FILE *fp = safe_fopen(name, "r", __FILE__, __LINE__);
    
    CHECK( fscanf(fp, "%f\n",  (real*   ) dz     ) );
    CHECK( fscanf(fp, "%f\n",  (real*   ) dx     ) );
    CHECK( fscanf(fp, "%f\n",  (real*   ) dy     ) );
    CHECK( fscanf(fp,  I"\n",  (integer*) dimmz  ) );
    CHECK( fscanf(fp,  I"\n",  (integer*) dimmx  ) );
    CHECK( fscanf(fp,  I"\n",  (integer*) dimmy  ) );
    CHECK( fscanf(fp, "%d\n",  (int*    ) nt_fwd ) );
    CHECK( fscanf(fp, "%d\n",  (int*    ) nt_bwd ) );
    CHECK( fscanf(fp, "%f\n",  (real*   ) dt     ) );
    CHECK( fscanf(fp, "%d\n",  (int*    ) stacki ) );
    
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
            log_error ( "Error while reading freqlist file"); 
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
            log_error ("Error while reading freqlist file"); 
            break;
        }   
        else if ( n == EOF )
        {
            break;
        }
    }
    fclose( freqfile ); 

    *nfreqs = count;

    log_info ( "There are %d wavelet frequencies to process)", *nfreqs );

    for( int i=0; i<count; i++)
        log_info ("     %.2f Hz", (*freqlist)[i] );
};

void* __malloc( size_t alignment, const integer size)
{
    void *buffer;
    int error;
    
    if( (error=posix_memalign( &buffer, alignment, size)) != 0) 
    {
        log_error ( "Cant allocate buffer correctly");
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
        log_error ( "%s:%d Cant open filename %s, openmode '%s'", srcfilename, linenumber, filename, mode);
        exit(-1);
    }   
    return temp;
};

void safe_fclose ( const char *filename, FILE* stream, char* srcfilename, int linenumber)
{
  if ( fclose( stream ) != 0)
  {
    log_error ( "%s:%d: Cant close file %s correctly!", srcfilename, linenumber, filename );
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
#ifndef DO_NOT_PERFORM_IO
    size_t res = fwrite( ptr, size, nmemb, stream);

    if( res != nmemb )
    {
        log_error ( "%s:%d: Error while fwrite", srcfilename, linenumber );
        abort();
    }
#endif
};

void safe_fread (void *ptr, size_t size, size_t nmemb, FILE *stream, char* srcfilename, int linenumber)
{
#ifndef DO_NOT_PERFORM_IO
    size_t res = fread( ptr, size, nmemb, stream);
    
    if( res != nmemb )
    {
        log_error ( "%s:%d: Error while fread", srcfilename, linenumber);
        abort();
    }
#endif
};

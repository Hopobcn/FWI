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

inline double dtime(void)
{
    double tseconds = 0.0;
    struct timeval mytime;
    gettimeofday( &mytime, (struct timezone*) 0);
    tseconds = (double) (mytime.tv_sec + (double) mytime.tv_usec * 1.0e-6);
    return (tseconds);
};

inline double TOGB(size_t bytes)
{
    return (bytes / (1024.f * 1024.f * 1024.f));
};

int max_int( int a, int b)
{
    return ((a >= b) ? a : b);
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


void read_fwi_parameters (const char *restrict fname,
                          real *lenz,
                          real *lenx,
                          real *leny,
                          real *vmin,
                          real *srclen,
                          real *rcvlen,
                          char *outputfolder)
{
    FILE *fp = safe_fopen(fname, "r", __FILE__, __LINE__ );

    IO_CHECK( fscanf( fp, "%f\n", (real*) lenz   ) );
    IO_CHECK( fscanf( fp, "%f\n", (real*) lenx   ) );
    IO_CHECK( fscanf( fp, "%f\n", (real*) leny   ) );
    IO_CHECK( fscanf( fp, "%f\n", (real*) vmin   ) );
    IO_CHECK( fscanf( fp, "%f\n", (real*) srclen ) );
    IO_CHECK( fscanf( fp, "%f\n", (real*) rcvlen ) );

    /* these three values are not needed for the shared memory implementation */
    int NotNeededValue;
    IO_CHECK( fscanf( fp, "%d\n", (int*) &NotNeededValue ) ); //nshots
    IO_CHECK( fscanf( fp, "%d\n", (int*) &NotNeededValue ) ); //ngrads
    IO_CHECK( fscanf( fp, "%d\n", (int*) &NotNeededValue ) ); //ntests
    IO_CHECK( fscanf( fp, "%d\n", (int*) &NotNeededValue ) ); //slavemem
    IO_CHECK( fscanf( fp, "%d\n", (int*) &NotNeededValue ) ); //workermem

    /* Recover the value of the output directory path */
    IO_CHECK( fscanf( fp, "%s\n",  outputfolder  ) );

    print_debug("Len (z,x,y) (%f,%f,%f) vmin %f scrlen %f rcvlen %f outputfolder '%s'",
      *lenz, *lenx, *leny, *vmin, *srclen, *rcvlen, outputfolder );

    fclose(fp);
};

void store_shot_parameters( const int     shotid,
                            const int     stacki,
                            const real    dt,
                            const int     nt_fwd,
                            const int     nt_bwd,
                            const real    dz,
                            const real    dx,
                            const real    dy,
                            const integer dimmz,
                            const integer dimmx,
                            const integer dimmy,
                            const char*   outputfolder,
                            const real    waveletFreq)
{
    char name[200];

    sprintf(name, "%s/shotparams_%2.1f.%05d.dat",
            outputfolder, waveletFreq, shotid);

    print_debug("Storing parameters for shot %d into %s", shotid, name);

    FILE *fp = safe_fopen(name, "w", __FILE__, __LINE__);

    fprintf(fp, "%f\n", dz     );
    fprintf(fp, "%f\n", dx     );
    fprintf(fp, "%f\n", dy     );
    fprintf(fp,  I"\n", dimmz  );
    fprintf(fp,  I"\n", dimmx  );
    fprintf(fp,  I"\n", dimmy  );
    fprintf(fp, "%d\n", nt_fwd );
    fprintf(fp, "%d\n", nt_bwd );
    fprintf(fp, "%f\n", dt     );
    fprintf(fp, "%d\n", stacki );

    fclose(fp);
};

void load_shot_parameters(const int   shotid,
                          int*        stacki,
                          real*       dt,
                          int*        nt_fwd,
                          int*        nt_bwd,
                          real*       dz,
                          real*       dx,
                          real*       dy,
                          integer*    dimmz,
                          integer*    dimmx,
                          integer*    dimmy,
                          const char* outputfolder,
                          const real  waveletFreq)
{
    char name[200];

    sprintf(name, "%s/shotparams_%2.1f.%05d.dat", outputfolder, waveletFreq, shotid);
    print_debug("Loading parameters for freq %d shot %d from %s", waveletFreq, shotid, name);

    FILE *fp = safe_fopen(name, "r", __FILE__, __LINE__);

    IO_CHECK( fscanf(fp, "%f\n",  (real*   ) dz     ) );
    IO_CHECK( fscanf(fp, "%f\n",  (real*   ) dx     ) );
    IO_CHECK( fscanf(fp, "%f\n",  (real*   ) dy     ) );
    IO_CHECK( fscanf(fp,  I"\n",  (integer*) dimmz  ) );
    IO_CHECK( fscanf(fp,  I"\n",  (integer*) dimmx  ) );
    IO_CHECK( fscanf(fp,  I"\n",  (integer*) dimmy  ) );
    IO_CHECK( fscanf(fp, "%d\n",  (int*    ) nt_fwd ) );
    IO_CHECK( fscanf(fp, "%d\n",  (int*    ) nt_bwd ) );
    IO_CHECK( fscanf(fp, "%f\n",  (real*   ) dt     ) );
    IO_CHECK( fscanf(fp, "%d\n",  (int*    ) stacki ) );

    safe_fclose( name, fp, __FILE__, __LINE__);
};


void load_freqlist( const char*  filename,
                          int*   nfreqs,
                          real** freqlist )
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
            print_error("Error while reading freqlist file");
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
            print_error("Error while reading freqlist file");
            break;
        }
        else if ( n == EOF )
        {
            break;
        }
    }
    fclose( freqfile );

    *nfreqs = count;

    print_info("A total of %d frequencies were found...", *nfreqs );
    for( int i=0; i<count; i++)
        print_info("     %.2f Hz", (*freqlist)[i] );
};


void create_output_volumes(const char*   outputfolder,
                           const integer VolumeMemory)
{
    print_debug("Creating output files in %s", outputfolder);

#ifdef DO_NOT_PERFORM_IO
    print_info("Warning: we are not doing any IO here.");
#else
    char fnamePrecond[300], fnameGradient[300];

    sprintf( fnameGradient, "%s/resultGradient.res", outputfolder);
    sprintf( fnamePrecond , "%s/resultPrecond.res", outputfolder);

    FILE *fGradient = safe_fopen( fnameGradient, "wb", __FILE__, __LINE__ );
    FILE *fPrecond  = safe_fopen( fnamePrecond , "wb", __FILE__, __LINE__ );

    int numIts = ceil( VolumeMemory / IO_CHUNK_SIZE );

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

void create_folder(const char* folder)
{
    if (mkdir_p(folder) != 0) {
        print_error("cant create folder %s (Err code: %s)", folder, strerror(errno));
        exit(-1);
    }
    print_debug("Folder '%s' created",folder);
};

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
                print_error("Error creating folder %s (Err code %s)", tmp, strerror(errno));
                return -1;
            }

            *p = '/';
        }
    }

    int rc = mkdir(tmp, S_IRWXU);
    if (rc != 0 && errno != EEXIST) {
        print_error("Error creating folder %s (Err code %s)", tmp, strerror(errno));
        return -1;
    }

    return 0;
}



FILE* safe_fopen(const char* filename,
                 const char* mode,
                 const char* srcfilename,
                 const int   linenumber)
{
    FILE* temp = fopen( filename, mode);

    if( temp == NULL)
    {
        print_error("Cant open filename %s, openmode '%s' (called from %s - %d)",
                    filename, mode, srcfilename, linenumber);
        exit(-1);
    }

    return temp;
};

void safe_fclose (const char* filename,
                        FILE* stream,
                  const char* srcfilename,
                  const int linenumber)
{
    if ( fclose( stream ) != 0)
    {
        print_error("Cant close filename %s (called from %s - %d)", filename, srcfilename, linenumber);
        abort();
    }
};


inline
void safe_fwrite (const void*  ptr,
                        size_t size,
                        size_t nmemb,
                        FILE*  stream,
                  const char*  srcfilename,
                  const int linenumber)
{
#ifdef DO_NOT_PERFORM_IO
    print_info("Warning: we are not doing any IO (called from %s).", __FUNCTION__);
#else
    if( stream == NULL )
    {
        print_error("Invalid stream\n");
        abort();
    }
    size_t res;

#if defined(LOG_IO_STATS)
    double start = dtime();
#endif
    res = fwrite( ptr, size, nmemb, stream);
#if defined(LOG_IO_STATS)
    double end = dtime() - start;

    double mbytes = (1.0 * size * nmemb) / (1024.0 * 1024.0);

    print_stats("Time %lf, elements %lu bytes %lu, MB %lf MB/s %lf",
                 end, nmemb, size*nmemb, mbytes, mbytes / end);
#endif

    if( res != nmemb )
    {
        print_error("Error while fwrite (called from %s - %d)", srcfilename, linenumber );
        abort();
    }
#endif
};

inline
void safe_fread (      void*  ptr,
                       size_t size,
                       size_t nmemb,
                       FILE*  stream,
                 const char*  srcfilename,
                 const int    linenumber)
{
#ifdef DO_NOT_PERFORM_IO
    print_info("Warning: we are not doing any IO (called from %s).", __FUNCTION__);
#else
    if( stream == NULL )
    {
        print_error("Invalid\n");
        abort();
    }

    size_t res = fread( ptr, size, nmemb, stream);

    if( res != nmemb )
    {
        print_error("Cant fread (called from %s - %d)", srcfilename, linenumber);
        print_error("Trying to read %lu elements, only %lu were recovered", nmemb, res);
        abort();
    }
#endif
};



void fwi_writelog(const char* SourceFileName,
                  const int   LineNumber,
                  const char* FunctionName,
                  const char* MessageHeader,
                  const char* fmt,
                  ...)
{
#if defined(USE_MPI)
    /* locate myself into the MPI world */
    int id;
    MPI_Comm_rank( MPI_COMM_WORLD, &id);
#else
    int id = 0;
#endif

    char LogFileName[50];
    sprintf(LogFileName, "fwi.%02d.log", id);

    FILE *fp = safe_fopen ( LogFileName, "a", __FILE__, __LINE__ );

    va_list args;
    va_start(args, fmt);
    fprintf(fp, "%s:%s:%d:%s: ", MessageHeader, SourceFileName, LineNumber, FunctionName );
    vfprintf(fp, fmt, args);
    fprintf(fp, "\n");
    va_end(args);

    safe_fclose ( LogFileName, fp, __FILE__, __LINE__);
};

void* __malloc( size_t alignment, const integer size)
{
    void *buffer;
    int error;

    if( (error=posix_memalign( &buffer, alignment, size)) != 0)
    {
        print_error("Cant allocate buffer correctly");
        abort();
    }

    return (buffer);
};

void __free ( void* ptr)
{
    free( ptr );
};


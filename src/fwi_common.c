#include "fwi_common.h"


/* extern variables declared in the header file */
const integer  WRITTEN_FIELDS =   12; /* >= 12.  */
const integer  HALO           =    4; /* >= 4    */ 
const integer  SIMD_LENGTH    =    8; /* # of real elements fitting into regs */
const real     IT_FACTOR      = 0.02;
const real     IO_CHUNK_SIZE  = 1024.f * 1024.f;

const size_t ALIGN_INT     = 16;
const size_t ALIGN_INTEGER = 16;
const size_t ALIGN_REAL    = 128; /* CUDA L1 cache line: 128bytes */

//extern 
FILE* logfile = NULL;



extent_t make_extent(size_t w, size_t h, size_t d)
{
    extent_t e = {0};

    e.width  = w;
    e.height = h;
    e.depth  = d;

    return e;
}


void log_info (const char *fmt, ...) 
{
#ifdef DEBUG

#if defined(USE_MPI)
    /* locate myself into the MPI world */
    int id;
    MPI_Comm_rank( MPI_COMM_WORLD, &id);
#else
    int id = 0;
#endif

    /* build log file name */
    char logname[50];
    sprintf( logname, "%02d.log", id);
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
    fprintf( flog, "%s: %s\n", timestr, str);
   
    /*  close file  */
    safe_fclose( logname, flog, __FILE__, __LINE__ );
#endif
};

void log_error (const char *fmt, ...) 
{ 
#if defined(USE_MPI)
    /* locate myself into the MPI world */
    int id;
    MPI_Comm_rank( MPI_COMM_WORLD, &id);
#else
    int id = 0;
#endif
    
    /* build log file name */
    char logname[50];
    sprintf( logname, "%02d.log", id);
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

void read_fwi_parameters (const char *fname,
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
    IO_CHECK( fscanf( fp, "%d\n", (int*) &NotNeededValue ) );
    IO_CHECK( fscanf( fp, "%d\n", (int*) &NotNeededValue ) );
    IO_CHECK( fscanf( fp, "%d\n", (int*) &NotNeededValue ) );
    IO_CHECK( fscanf( fp, "%d\n", (int*) &NotNeededValue ) );
    IO_CHECK( fscanf( fp, "%d\n", (int*) &NotNeededValue ) );
 
    /* Recover the value of the output directory path */
    IO_CHECK( fscanf( fp, "%s\n",  outputfolder  ) );

    print_debug("Len (z,x,y) (%f,%f,%f) vmin %f scrlen %f rcvlen %f outputfolder '%s'",
      *lenz, *lenx, *leny, *vmin, *srclen, *rcvlen, outputfolder );

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
        print_error("cant create folder %s (Err code: %s)", folder, strerror(errno));
        exit(-1);
    }
    print_debug("Folder '%s' created",folder);
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

void store_shot_parameters(int     shotid,
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
                           char    *outputfolder, 
                           real    waveletFreq)
{
    char name[200];

    sprintf(name, "%s/shotparams_%2.1f.%05d.dat", 
            outputfolder, waveletFreq, shotid);

    print_debug("Storing parameters for shot %d into %s", shotid, name);

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

void load_shot_parameters(int     shotid,
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
                          char    *outputfolder,
                          real    waveletFreq)
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

void* malloc3d_host(dim_t* dim, const size_t alignment, extent_t req)
{
    const size_t ALIGN = alignment;

    void* h_base_ptr;

    size_t padding = ( (size_t)abs((int)ALIGN - (int)(req.width*sizeof(real))) )%ALIGN;
    assert( padding%sizeof(real) == 0 );
    
    dim->pitch = req.width + (padding / sizeof(real));
    dim->zsize = req.width;
    dim->xsize = req.height;
    dim->ysize = req.depth;

    const size_t size = (dim->pitch * dim->xsize * dim->ysize + HALO) * sizeof(real);
    /* 
     * 'size' is a perfect aligned 3D volume size [pitch x xsize x ysize] if your access have stride 0
     *  BUT our accesses have a stride HALO: ptr[HALO] is typically the first internal position.
     *      so we have to align to that stride also!
     */

    h_base_ptr = (void*) malloc( size+ALIGN-1 );

    uintptr_t mask = ~(uintptr_t)(ALIGN-1);

    void** h_ptr = (void**) ((((uintptr_t)h_base_ptr+ALIGN-1) & mask) + (ALIGN-4-12));

#if defined(DEBUG)
    fprintf(stdout, "--malloc3d: h_base_ptr %p h_alig_ptr %p\n", 
            h_base_ptr, *h_ptr);
#endif

    /* 
     * since we have to call 'free' and 'acc_free' with base pointers 
     * we have to store them in the innaccessible part 
     * WARNING: any ovewrite of this information will cause a SEGFAULT in the free
     */
    h_ptr[-1] = h_base_ptr;

    return h_ptr;
}

void free3d_host(void** h_align_ptr)
{
    /* deallocate memory from the base pointers */
    void* h_base_ptr = h_align_ptr[-1];

    free(h_base_ptr);
}

#if defined(_OPENACC)
void* malloc3d_device(dim_t* dim, const size_t alignment, extent_t req, void** h_align_ptr)
{
    const size_t ALIGN = alignment;

    void* d_base_ptr;

    size_t padding = ( (size_t)abs((int)ALIGN - (int)(req.width*sizeof(real))) )%ALIGN;
    assert( padding%sizeof(real) == 0 );
    
    dim->pitch = req.width + (padding / sizeof(real));
    dim->zsize = req.width;
    dim->xsize = req.height;
    dim->ysize = req.depth;

    const size_t size = (dim->pitch * dim->xsize * dim->ysize + HALO) * sizeof(real);
    /* 
     * 'size' is a perfect aligned 3D volume size [pitch x xsize x ysize] if your access have stride 0
     *  BUT our accesses have a stride HALO: ptr[HALO] is typically the first internal position.
     *      so we have to align to that stride also!
     */

    d_base_ptr =     acc_malloc( size+ALIGN-1 );

    uintptr_t mask = ~(uintptr_t)(ALIGN-1);

    // TODO: figure out what this -4-12 is..
    void** d_ptr = (void**) ((((uintptr_t)d_base_ptr+ALIGN-1) & mask) + (ALIGN-4-12));

#if defined(DEBUG)
    fprintf(stdout, "--malloc3d: d_base_ptr %p d_alig_ptr %p\n", 
            d_base_ptr, *d_ptr);
#endif

    /* 
     * since we have to call 'free' and 'acc_free' with base pointers 
     * we have to store them in the innaccessible part 
     * WARNING: any ovewrite of this information will cause a SEGFAULT in the free
     */
    h_align_ptr[-2] = d_base_ptr;

    return d_ptr;
}

void free3d_device(void** h_align_ptr)
{
    /* deallocate memory from the base pointers */
    void* d_base_ptr = h_align_ptr[-2];

    acc_free(d_base_ptr);
}
#endif


FILE* safe_fopen(const char *filename, const char *mode, const char* srcfilename, const int linenumber)
{
    FILE* temp = fopen( filename, mode);
    
    if( temp == NULL){
        print_error("Cant open filename %s, openmode '%s' (called from %s - %d)", 
                    filename, mode, srcfilename, linenumber);
        exit(-1);
    }
    return temp;
};

void safe_fclose ( const char *filename, FILE* stream, const char* srcfilename, const int linenumber)
{
    if ( fclose( stream ) != 0)
    {
        print_error("Cant close filename %s (called from %s - %d)", filename, srcfilename, linenumber);
        abort();
    }

/*if ( unlink(filename)  != 0)
  {
    fprintf(stderr, "%s:%d: Cant unlink file %s correctly!\n", srcfilename, linenumber, filename );
    abort();
  }*/
};


inline 
void safe_fwrite (const void *ptr, size_t size, size_t nmemb, FILE *stream, const char* srcfilename, const int linenumber)
{
#ifdef DO_NOT_PERFORM_IO
    print_info("Warning: we are not doing any IO (called from %s).", __FUNCTION__);
#else
    if( stream == NULL ){
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
void safe_fread (void *ptr, size_t size, size_t nmemb, FILE *stream, const char* srcfilename, const int linenumber)
{
#ifdef DO_NOT_PERFORM_IO
    print_info("Warning: we are not doing any IO (called from %s).", __FUNCTION__);
#else
    if( stream == NULL ){
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



void fwi_writelog(const char *SourceFileName, 
                  const int LineNumber,
                  const char *FunctionName,
                  const char* MessageHeader,
                  const char *fmt,
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


int parse_env(const char* name)
{
    char* value = getenv(name);
    if (value != NULL)
    {
        return atoi(value);
    }
    return 0;
}

#if defined(USE_MPI)
int mpi_get_rank()
{
#if defined(OPEN_MPI)
    int rank       = parse_env("OMPI_COMM_WORLD_RANK");
#elif defined(MPICH)
    int rank       = parse_env("MV2_COMM_WORLD_RANK");
#else
    int rank       = 0;
#endif
    return rank;
}

int mpi_get_local_rank()
{
#if defined(OPEN_MPI)
    int local_rank = parse_env("OMPI_COMM_WORLD_LOCAL_RANK");
#elif defined(MPICH)
    int local_rank = parse_env("MV2_COMM_WORLD_LOCAL_RANK");
#else
    int local_rank = 0;
#endif
}
#endif

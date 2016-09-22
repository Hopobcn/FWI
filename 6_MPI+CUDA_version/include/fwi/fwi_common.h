/*
 * =====================================================================================
 *
 *       Filename:  fwi_common.h
 *
 *    Description:
 *
 *        Version:  1.0
 *        Created:  10/12/15 10:34:04
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  YOUR NAME (),
 *   Organization:
 *
 * =====================================================================================
 */

#ifndef _FWI_COMMON_H_
#define _FWI_COMMON_H_

// When included before <stdlib.h>, solves implicit declaration of posix_memalign()
// http://stackoverflow.com/questions/32438554/warning-implicit-declaration-of-posix-memalign
#define _POSIX_C_SOURCE 200809L
#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>
#include <stdint.h>
#include <assert.h>
#include <math.h>
#include <sys/time.h>
#include <string.h>
#include <errno.h>
#include <sys/stat.h>
#include <unistd.h>
#include <time.h>

#if defined(USE_MPI)
#include <mpi.h>
#endif

#if defined(_OPENMP)
#include <omp.h>
#endif

#if defined(_OPENACC) && !defined(__NVCC__)
#include <openacc.h>
#endif

#if defined(TRACE_CUDA)
#include <nvToolsExt.h>
#endif


/* data types definition */
typedef float  real;
typedef int    integer;

#define I "%d"     // integer printf symbol

typedef enum {RTM_KERNEL, FM_KERNEL} propagator_t;
typedef enum {FORWARD   , BACKWARD, FWMODEL}  time_d;

/* simulation parameters */
extern const integer WRITTEN_FIELDS;
extern const integer HALO;
extern const integer SIMD_LENGTH;
extern const real    IT_FACTOR;
extern const real    IO_CHUNK_SIZE;

extern const size_t ALIGN_INT;
extern const size_t ALIGN_INTEGER;
extern const size_t ALIGN_REAL;

/* (MPI-local) file for logging */
extern FILE* logfile;

double TOGB(size_t bytes);

/*  Compiler compatiblity macros */
#if defined(__GNUC__)
  /* http://stackoverflow.com/questions/25667901/assume-clause-in-gcc*/ \
    #define __assume(_cond) do { if (!(_cond)) __builtin_unreachable(); } while (0)
#endif

/*  Compiler macro to suppress unused variable warnings */
#if defined(UNUSED)
#elif defined(__GNUC__)
    #define UNUSED(x) (x) __attribute__((unused))
#else
    #define UNUSED(x) x
#endif

#define CHECK(error) { checkErrors((error), __FILE__, __LINE__); }
static inline void checkErrors(const integer error, const char *filename, int line)
{
    if ( error < 0 ) {                     
        fprintf(stderr, "ERROR: %d in %s:%d\n", error, filename, line);
        exit(-1);
    }
};

FILE* safe_fopen  ( const char *filename, const char *mode, const char* srcfilename, const int linenumber);
void  safe_fclose ( const char *filename, FILE* stream, const char* srcfilename, const int linenumber);
void  safe_fwrite ( const void *ptr, size_t size, size_t nmemb, FILE *stream, const char* srcfilename, const int linenumber );
void  safe_fread  (       void *ptr, size_t size, size_t nmemb, FILE *stream, const char* srcfilename, const int linenumber );
integer roundup(integer number, integer multiple);

void log_info  (const char *fmt, ...);
void log_error (const char *fmt, ...);

int max_int( int a, int b);

double dtime(void);

void read_fwi_parameters (const char *fname,
                          real *lenz,
                          real *lenx,
                          real *leny,
                          real *vmin,
                          real *srclen,
                          real *rcvlen,
                          char *outputfolder);

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
                           real    waveletFreq);

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
                          real    waveletFreq);

void load_freqlist (  const char*  filename,
                            int*   nfreqs,
                            real** freqlist );

void* __malloc ( const size_t alignment, const integer size);
void  __free   ( void *ptr );

void create_output_volumes(char* outputfolder, integer VolumeMemory);

int mkdir_p(const char *dir);

void create_folder(const char *folder);


#define print_error(M, ...)     fwi_writelog(__FILE__, __LINE__, __func__, "ERROR ", M, ##__VA_ARGS__)
#define print_info(M, ...)      fwi_writelog(__FILE__, __LINE__, __func__, "INFO  ", M, ##__VA_ARGS__)
#define print_stats(M, ...)     fwi_writelog(__FILE__, __LINE__, __func__, "STATS ", M, ##__VA_ARGS__)

#if defined(DEBUG)
  #define print_debug(M, ...)  fwi_writelog(__FILE__, __LINE__, __func__, "DEBUG ", M, ##__VA_ARGS__)
#else
  #define print_debug(M, ...)
#endif

void fwi_writelog(const char *SourceFileName, 
                  const int LineNumber,
                  const char *FunctionName,
                  const char* MessageHeader,
                  const char *fmt,
                  ...);

#if defined(TRACE_CUDA)
    #define PUSH_RANGE nvtxRangePush(__func__);
    #define POP_RANGE  nvtxRangePop();
#else
    #define PUSH_RANGE
    #define POP_RANGE
#endif

int parse_env(const char* name);

// 
// GPU-Affinity related functions:
//
#if defined(USE_MPI)
int mpi_get_rank();
int mpi_get_local_rank();
#endif
//Obs: due a BUG with PGI 16.5 & CUDA 7.5 we can't just 'include <cuda.h>' in C files
//     we have to use NVCC and PGI separately and link the result
#if 0
#if defined(__cplusplus)
extern "C"
#endif
int select_gpu_and_pin_proc(int rank, int local_rank);
#endif

#endif // end of _FWI_COMMON_H_ definition

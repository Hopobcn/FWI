/*
 * =============================================================================
 * Copyright (c) 2016, Barcelona Supercomputing Center (BSC)
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *     * Redistributions of source code must retain the above copyright
 *       notice, this list of conditions and the following disclaimer.
 *     * Redistributions in binary form must reproduce the above copyright
 *       notice, this list of conditions and the following disclaimer in the
 *       documentation and/or other materials provided with the distribution.
 *     * Neither the name of the <organization> nor the
 *       names of its contributors may be used to endorse or promote products
 *       derived from this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL <COPYRIGHT HOLDER> BE LIABLE FOR ANY
 * DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 * =============================================================================
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

#define IO_CHECK(error) { checkErrors((error), __FILE__, __LINE__); }
static inline void checkErrors(const integer error, const char *filename, int line)
{
    if ( error < 0 ) {
        fprintf(stderr, "ERROR: %d in %s:%d\n", error, filename, line);
        exit(-1);
    }
};

/*
 * \brief Computes the number of seconds since the Epoch used by gettimeofday
 * \return Seconds since Epoch with microsecond resolution
 */
double dtime(void);

/*
 * \brief Convert from Bytes to Giga-Bytes
 * \param bytes [in] number of bytes
 * \return number of GB
 */
double TOGB(size_t bytes);

/*
 * \brief Computes the maximum value between two integers
 * \param a [in] integer
 * \param b [in] integer
 * \return maximum integer between a and b
 */
int max_int( int a, int b);

/*
 * This function is intended to round up a number (number) to the nearest multiple of the register
 * size. In this way, we assure that the dimensions of the domain are suited to the most aggressive
 * compiler optimizations.
 *
 * \param number   [in] number to be roundup
 * \param multiple [in] next
 */
integer roundup(integer number, integer multiple);


/*
 * \brief reads FWI parameters.
 *
 *  Opens a filename file and reads all parameters but only some are used in FWI.
 *  For mor information of the parameters layout go to data/README.md
 *  If debug is enabled prints those parameters
 *  in a log file.
 *
 * \param fname        [in] filename
 * \param lenz         [out] length in Z dimension (meters)
 * \param lenx         [out] length in X dimension (meters)
 * \param leny         [out] length in Y dimension (meters)
 * \param vmin         [out] velocity (m/s) (used to compute dz/dx/dy)
 * \param srclen       [out] (used to compute forw_steps )
 * \param rcvlen       [out] (used to compute back_steps)
 * \param outputfolder [out] name of the folder were IO files is put
 */
void read_fwi_parameters (const char *restrict fname,
                          real *lenz,
                          real *lenx,
                          real *leny,
                          real *vmin,
                          real *srclen,
                          real *rcvlen,
                          char *outputfolder);

/*
 * \brief stores shot parameters for later IO post-processing
 *
 *  Opens a filename file which name is parametrized by:
 *    - outputfolder
 *    - waveletFreq
 *    - shotid
 *  so we can identify different IO files by their freq & shotid
 *
 * \param shotid       [in]
 * \param stacki       [in]
 * \param dt           [in]
 * \param nt_fwd       [in]
 * \param nt_bwd       [in]
 * \param dz           [in]
 * \param dx           [in]
 * \param dy           [in]
 * \param dimmz        [in]
 * \param dimmx        [in]
 * \param dimmy        [in]
 * \param outputfolder [in]
 * \param waveletFreq  [in]
 */
void store_shot_parameters (const int     shotid,
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
                            const real    waveletFreq);

/*
 * \brief load shot parameters previously stored with 'store_shot_parameters'
 *
 *  Opens a filename file which name is parametrized by:
 *    - outputfolder
 *    - waveletFreq
 *    - shotid
 *  and initializes all parameters
 *
 * \param shotid       [in]
 * \param stacki       [in]
 * \param dt           [in]
 * \param nt_fwd       [in]
 * \param nt_bwd       [in]
 * \param dz           [in]
 * \param dx           [in]
 * \param dy           [in]
 * \param dimmz        [in]
 * \param dimmx        [in]
 * \param dimmy        [in]
 * \param outputfolder [in]
 * \param waveletFreq  [in]
 */
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
                          const real  waveletFreq);

/*
 * \brief load shot parameters previously stored with 'store_shot_parameters'
 *
 * \param filename [in]  frequency filename
 * \param nfreqs   [out] number of frequencies
 * \param freqlist [out] list of frequencies
 */
void load_freqlist (const char*  filename,
                          int*   nfreqs,
                          real** freqlist );

/*
 * \brief creates files to store 3D volumes contents (IO)
 *
 * Creates files to store the final preconditioner and gradient results
 *
 * \param outputfolder [in]  folder where snapshot data is stored
 * \param VolumeMemory [in]  memory needed to store the domain
 */
void create_output_volumes(const char*   outputfolder,
                           const integer VolumeMemory);

/*
 * \brief creates a folder
 *
 * During execution FWI creates temporal folders to organize necessary data
 * for the execution
 *
 * \param folder [in] name of the folder to create
 */
void create_folder(const char *folder);

/*
 * \brief creates the hierarchy of folders requested, if they do not exist
 *
 * \param dir [in] name of the directory to create
 * \return 0 if successful, 0 otherwise
 */
int mkdir_p(const char *dir);

/*
 * \brief opens a file
 *
 * \param filename    [in] name of the file to be opened
 * \param mode        [in] file access mode (r/w/a/r+/w+/a+)
 * \param srcfilename [in] source file name
 * \param linenumber  [in] source file line number
 * \return file descriptor of that filename
 */
FILE* safe_fopen (const char* filename,
                  const char* mode,
                  const char* srcfilename,
                  const int   linenumber);

/*
 * \brief closes a file
 *
 * \param filename    [in] name of the file to be closed
 * \param stream      [in] stream to be closed
 * \param srcfilename [in] source file name
 * \param linenumber  [in] source file line number
 */
void safe_fclose (const char* filename,
                        FILE* stream,
                  const char* srcfilename,
                  const int   linenumber);

/*
 * \brief write a binary buffer to file
 *
 * \param ptr         [in] base pointer of an array to be written
 * \param size        [in] size of each element
 * \param nmemb       [in] number of elements to be written
 * \param srcfilename [in] source file name
 * \param linenumber  [in] source file line number
 */
void safe_fwrite (const void*  ptr,
                        size_t size,
                        size_t nmemb,
                        FILE*  stream,
                  const char*  srcfilename,
                  const int    linenumber);

/*
 * \brief reads a binary buffer from file and store its conters in memory
 *
 * \param ptr         [in] base pointer of an array to be read
 * \param size        [in] size of each element
 * \param nmemb       [in] number of elements to be written
 * \param srcfilename [in] source file name
 * \param linenumber  [in] source file line number
 */
void safe_fread  (      void*  ptr,
                        size_t size,
                        size_t nmemb,
                        FILE*  stream,
                  const char*  srcfilename,
                  const int    linenumber);


#define print_error(M, ...)     fwi_writelog(__FILE__, __LINE__, __func__, "ERROR ", M, ##__VA_ARGS__)
#define print_info(M, ...)      fwi_writelog(__FILE__, __LINE__, __func__, "INFO  ", M, ##__VA_ARGS__)
#define print_stats(M, ...)     fwi_writelog(__FILE__, __LINE__, __func__, "STATS ", M, ##__VA_ARGS__)

#if defined(DEBUG)
  #define print_debug(M, ...)  fwi_writelog(__FILE__, __LINE__, __func__, "DEBUG ", M, ##__VA_ARGS__)
#else
  #define print_debug(M, ...)
#endif

/*
 * \brief C variadic function that writes a variable number of arguments to a file
 *
 * \param SourceFileName [in] source file name
 * \param LineNumber     [in] source file line number
 * \param FunctionName   [in] function name
 * \param MessageHeader  [in] message handler (error/info/stats/debug)
 * \param fm             [in] argument list
 * \param ...            [in] ellipsis
 */
void fwi_writelog(const char* SourceFileName,
                  const int   LineNumber,
                  const char* FunctionName,
                  const char* MessageHeader,
                  const char* fmt,
                  ...);

#if defined(TRACE_CUDA)
    #define PUSH_RANGE nvtxRangePush(__func__);
    #define POP_RANGE  nvtxRangePop();
#else
    #define PUSH_RANGE
    #define POP_RANGE
#endif

/*
 * \brief allocate aligned memory
 *
 * \param alignment [in]
 * \param size      [in]
 * \return base pointer of the allocated region
 */
void* __malloc ( const size_t alignment, const integer size);

/*
 * \brief free memory
 *
 * \param ptr [in] base pointer of the allocation to be freed
 */
void  __free   ( void *ptr );

#endif // end of _FWI_COMMON_H_ definition

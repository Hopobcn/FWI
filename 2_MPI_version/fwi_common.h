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
#include <math.h>
#include <sys/time.h>
#include <string.h>
#include <errno.h>
#include <sys/stat.h>
#include <unistd.h>
#include <time.h>
#include <mpi.h>
#include <omp.h>


/* data types definition */
typedef float  real;
typedef int    integer;

#define I "%d"     // integer printf symbol

typedef enum {RTM_KERNEL, FM_KERNEL} propagator_t;
typedef enum {FORWARD   , BACKWARD, FWMODEL}  time_d;

/* simulation parameters */
extern const integer WRITTEN_FIELDS;
extern const integer HALO;
extern const real    IT_FACTOR;
extern const real    IO_CHUNK_SIZE;

extern const size_t ALIGN_INT;
extern const size_t ALIGN_INTEGER;
extern const size_t ALIGN_REAL;

/* (MPI-local) file for logging */
extern FILE* logfile;


#define TOGB(bytes) (bytes / (1024.f * 1024.f * 1024.f))


FILE* safe_fopen  ( const char *filename, char *mode, char* srcfilename, int linenumber);
void  safe_fclose ( const char *filename, FILE* stream, char* srcfilename, int linenumber);
void  safe_fwrite ( void *ptr, size_t size, size_t nmemb, FILE *stream, char* srcfilename, int linenumber );
void  safe_fread  ( void *ptr, size_t size, size_t nmemb, FILE *stream, char* srcfilename, int linenumber );

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

void store_shot_parameters( int     shotid,
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
                           char    *outputfolder);

void load_shot_parameters( int     shotid,
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
                          char    *outputfolder);

void load_freqlist (  const char* filename, 
											int *nfreqs, 
											real **freqlist ); 

void* __malloc ( const size_t alignment, const integer size);
void  __free   ( void *ptr );

void create_output_volumes(char* outputfolder, integer VolumeMemory);

int mkdir_p(const char *dir); 
void create_folder(const char *folder);

#endif // end of _FWI_COMMON_H_ definition


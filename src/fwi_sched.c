/*
 * =============================================================================
 * Copyright (c) 2016-2018, Barcelona Supercomputing Center (BSC)
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

#include "fwi/fwi_sched.h"

void schedule_free( schedule_t s )
{
    free ( s.freq );
    free ( s.forws );
    free ( s.backs );
    free ( s.stacki );
    free ( s.dt );
    free ( s.dz );
    free ( s.dx );
    free ( s.dy );
    free ( s.dimmz );
    free ( s.dimmx );
    free ( s.dimmy );
    free ( s.ppd );
    free ( s.nworkers );
};

schedule_t load_schedule( const char* filename ) 
{
    /* open file to store the schedule */
    char* fwipath = read_env_variable("FWIDIR");
    int pathlen = strlen(fwipath) + 200;
    char filepath[pathlen];
    sprintf(filepath, "%s/data/%s", fwipath, filename);

    FILE* fschedule = safe_fopen( filepath, "r", __FILE__, __LINE__);

    if ( fschedule == NULL ){
        fprintf(stderr, "Cant open scheduling file %s\n", filepath );
        abort();
    }

    schedule_t s;

    /* read the number of frequencies to be processed */
    if (fscanf( fschedule, "%d\n", &s.nfreqs ) != 1)
        print_error("cannot read nfreqs from file");
    fprintf(stderr, "Number of frequencies %d\n", s.nfreqs);

    /* read the number of shots to be processed */
    if (fscanf( fschedule, "%d\n", &s.nshots ) != 1)
        print_error("cannot read nshots from file");
    fprintf(stderr, "Number of shots %d\n", s.nshots);

    /* read the number of gradient iterations */
    if (fscanf( fschedule, "%d\n", &s.ngrads ) != 1)
        print_error("cannot read ngrads from file");
    fprintf(stderr, "Number of gradient iterations %d\n", s.ngrads);

    /* read the number of test iterations */
    if (fscanf( fschedule, "%d\n", &s.ntests ) != 1)
        print_error("cannot read ntests from file");
    fprintf(stderr, "Number of test iterations %d\n", s.ntests);

    /* read the name of the output directory */
    if (fscanf( fschedule, "%s\n", s.outputfolder ) != 1)
        print_error("cannot read outputfolder from file");
    fprintf(stderr, "Output directory path: %s\n", s.outputfolder );

    /* allocate memory for the rest of the parameters */
    s.freq     = (real*   ) malloc( s.nfreqs * sizeof(real   ));
    s.forws    = (integer*) malloc( s.nfreqs * sizeof(integer));
    s.backs    = (integer*) malloc( s.nfreqs * sizeof(integer));
    s.stacki   = (integer*) malloc( s.nfreqs * sizeof(integer));
    s.dt       = (real*   ) malloc( s.nfreqs * sizeof(real   ));
    s.dz       = (real*   ) malloc( s.nfreqs * sizeof(real   ));
    s.dx       = (real*   ) malloc( s.nfreqs * sizeof(real   ));
    s.dy       = (real*   ) malloc( s.nfreqs * sizeof(real   ));
    s.dimmz    = (integer*) malloc( s.nfreqs * sizeof(integer));
    s.dimmx    = (integer*) malloc( s.nfreqs * sizeof(integer));
    s.dimmy    = (integer*) malloc( s.nfreqs * sizeof(integer));
    s.ppd      = (integer*) malloc( s.nfreqs * sizeof(integer));
    s.nworkers = (integer*) malloc( s.nfreqs * sizeof(integer));

    /* read the list of simulation parameters */
    for( int i=0; i < s.nfreqs; i++ ) 
    {
        if (fscanf( fschedule, "%f %d %d %d %f %f %f %f %d %d %d %d %d\n",
                &s.freq[i], &s.forws[i], &s.backs[i], &s.stacki[i],
                &s.dt[i], &s.dz[i], &s.dy[i], &s.dx[i],
                &s.dimmz[i], &s.dimmx[i], &s.dimmy[i],
                &s.ppd[i], &s.nworkers[i]) != 13)
            print_error("cannot read simulation parameters from file");

#if defined(SHARED_MEMORY_RUN)
        if ( s.nworkers[i] != 1 ) {
            print_error("This execution requieres more than one worker to be computed. Please check the schedule file.");
            abort();
        }
#endif
    }

    /* show up a summary */
    print_info("Info: number of frequencies %d, number of shots %d", s.nfreqs, s.nshots);

    for( int i=0; i < s.nfreqs; i++ )
    {
        print_info("%f %d %d %d %f %f %f %f %d %d %d %d %d", 
            s.freq[i], s.forws[i], s.backs[i], s.stacki[i], 
            s.dt[i], s.dz[i], s.dy[i], s.dx[i], 
            s.dimmz[i], s.dimmx[i], s.dimmy[i], 
            s.ppd[i], s.nworkers[i]);
    }
    /* clean up and resume */
    fclose( fschedule ); fschedule = NULL;

    return (s);
};


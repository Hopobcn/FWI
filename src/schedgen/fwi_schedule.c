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

/*
 * Estimates the number of workers (accelators) needed to carry on the
 * simulation. To that end, we use the max. wavelet frequency of the
 * simulation and the available memory on the accelerator device (in GB).
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <math.h>

#include "fwi/fwi_sched.h"

const double BytesInGB = 1024.f * 1024.f * 1024.f;

int main(int argc, char *argv[])
{
    /*
     * Check input arguments.
     */
    if ( argc != 3 ) {
        printf("Usage: %s <paramfile> <freqsfile>", argv[0]);
        print_error("Invalid argument count. Parameters and frequencies files are requiered.");
        abort();
    }

    if ( argv[1] == NULL || argv[2] == NULL )
    {
        print_error("Invalid arguments");
        abort();
    }

    /*
     * Local parameters that define the simulation
     */
    int nshots;
    int ngrads;
    int ntests;
    int nfreqs;
    real slavemem;
    real workermem;
    int nworkers;
    real lenz;
    real lenx;
    real leny;
    real vmin;
    real srclen;
    real rcvlen;
    real *frequencies;
    size_t device_mem;
    char outputfolder[1024];
    char originalfolder[1024];

    /*
     * Modify input parameters according to FWIDIR value
     */
    char* fwipath = read_env_variable("FWIDIR");
    int pathlen = strlen(fwipath) + 200;

    char paramfile[pathlen];
    char freqsfile[pathlen];
    char schedfile[pathlen];

    sprintf( paramfile,  "%s/data/%s", fwipath, argv[1]);
    sprintf( freqsfile,  "%s/data/%s", fwipath, argv[2]);
    sprintf( schedfile,  "%s/data/fwi_schedule.txt", fwipath );

    read_fwi_parameters ( paramfile,
        &lenz,
        &lenx,
        &leny,
        &vmin,
        &srclen,
        &rcvlen,
        &nshots,
        &ngrads,
        &ntests,
        &slavemem,
        &workermem,
        originalfolder);

    printf("Path provided in fwi_params.txt is: %s\n", originalfolder);

    IO_CHECK(sprintf(outputfolder, "%s", originalfolder));
    printf("Using default global parallel file system: %s\n", outputfolder);

    load_freqlist( freqsfile, &nfreqs, &frequencies );

    /* open file to store the schedule */
    FILE* fschedule = safe_fopen( schedfile, "w", __FILE__, __LINE__);

    if ( fschedule == NULL ) {
        print_error("ERROR: Cant create scheduling file %s", schedfile);
        abort();
    }

    /* store the number of frequencies and shots on the schedule file */
    IO_CHECK( fprintf( fschedule, "%d\n", nfreqs) );
    IO_CHECK( fprintf( fschedule, "%d\n", nshots) );
    IO_CHECK( fprintf( fschedule, "%d\n", ngrads) );
    IO_CHECK( fprintf( fschedule, "%d\n", ntests) );
    IO_CHECK( fprintf( fschedule, "%s\n", outputfolder ));

    for( int freq = 0; freq < nfreqs; freq++)
    {
        /* our calculations are based on the max frequency */
        real waveletFreq = frequencies[freq];
        print_info("Estimating resources for %f Hz...", waveletFreq);

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

        print_info("Value of stacki %d forward steps %d backward steps %d dt %f", 
                    stacki, forw_steps, back_steps, dt );
        print_info("Value of dz = %f dx = %f dy = %f", dz, dx, dy);

        integer yplane_mem = dimmx * dimmz * 58 * sizeof(real); 


        /* if the accelerator memory is 0 means that it wont be used */
        if ( workermem ) { device_mem = (0.8 * workermem) * BytesInGB; }
        else             { device_mem = (0.8 * slavemem ) * BytesInGB; }

        /* compute the number of y-planes fitting into the device */
        integer ppd = device_mem / yplane_mem;

        /* minum number of planes have to be satisfied */
        if ( ppd < 4 * HALO ){
            print_error("ERROR: At least 4*HALO planes must fit into a single node.");
            print_error("ERROR: You wont be able to run this configuration..");
        }

        /* compute how many nodes are required */
        if ( dimmy < ppd ){
            nworkers = 1;
        } else {
            nworkers = dimmy / ppd;

            if ( dimmy % ppd != 0)
                nworkers += 1;
        }

        printf("  ------------------------------------------------------------------------------------\n\n");
        printf("                                FWI SCHEDULE GENERATOR\n");
        printf("	There are %d y-planes to compute, each worker can hold %d of these planes\n", dimmy, ppd);
        printf("	At this frequency (%f Hz) we'll need %d workers per shot\n", waveletFreq, nworkers );
        printf("	There are %d shots to be computed, so %d slave nodes are needed.\n", nshots, nshots );
        printf("  ------------------------------------------------------------------------------------\n\n");

        IO_CHECK( fprintf( fschedule, "%f %d %d %d %f %f %f %f %d %d %d %d %d\n", 
        waveletFreq, forw_steps, back_steps, stacki, dt, dz, dy, dx, dimmz, dimmx, dimmy, ppd, nworkers));
    }

    fclose( fschedule ); fschedule = NULL;

    print_info("End of model generator routine! STATUS: OK");

    return 0;
};


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

#include "fwi/fwi_kernel.h"


int main(int argc, const char *argv[])
{
    if (argc < 3) {
        fprintf(stderr, "Usage: %s <params_file> <frequency_file>\n", argv[0]);
        exit(0);
    }

    /* set seed for random number generator */
    srand(314);

    real lenz,lenx,leny,vmin,srclen,rcvlen;
    char outputfolder[200];

    fprintf(stderr, "Loading parameter from %s file\n", argv[1]);
    read_fwi_parameters( argv[1], &lenz, &lenx, &leny, &vmin, &srclen, &rcvlen, outputfolder);

    /* create synthetic velocity model */
    int nfreqs;
    real *frequencies;

    load_freqlist( argv[2], &nfreqs, &frequencies);

    for(int i=0; i<nfreqs; i++)
    {
        real waveletFreq = frequencies[i];
        fprintf(stderr, "Creating synthetic velocity input model for %f Hz freq\n", waveletFreq );

        /* compute discretization deltas, 16 == puntos por longitud de onda */
        real dx = vmin / (16.0 * waveletFreq);
        real dy = vmin / (16.0 * waveletFreq);
        real dz = vmin / (16.0 * waveletFreq);

        /* number of cells along axis */
        integer dimmz = roundup( ceil( lenz / dz ) + 2*HALO, HALO);
        integer dimmy = roundup( ceil( leny / dy ) + 2*HALO, HALO);
        integer dimmx = roundup( ceil( lenx / dx ) + 2*HALO, HALO);

        const integer numberOfCells = dimmz * dimmy * dimmx;

        fprintf(stderr, "Elements/array = "I"\n", numberOfCells);

        char modelname[300];
        sprintf( modelname, "../data/inputmodels/velocitymodel_%.2f.bin", waveletFreq );

        FILE* model = safe_fopen( modelname, "wb", __FILE__, __LINE__);

        real *buffer = __malloc( ALIGN_REAL, sizeof(real) * numberOfCells);

        /* safe dummy buffer */
        for(int i = 0; i < WRITTEN_FIELDS; i++)
        {
            /* fill the buffer with random numbers in [-1,1] interval */
            //for(int j = 0; j < numberOfCells; j++)
            //    buffer[j] = (rand() % 2) -1.0;     <<--- provoca NaNs en l'execucio
            set_array_to_random_real( buffer, numberOfCells );

            safe_fwrite( buffer, sizeof(real), numberOfCells, model, __FILE__, __LINE__);
        }

        /* free buffer */
        __free( buffer );

        /*  close model file */
        safe_fclose( modelname, model, __FILE__, __LINE__);

        fprintf(stderr, "Model %s created correctly\n", modelname);
    }

    __free( frequencies );
    fprintf(stderr, "End of the program\n");
    return 0;
}

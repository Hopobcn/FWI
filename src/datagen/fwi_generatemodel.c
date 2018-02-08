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
#include "fwi/fwi_kernel.h"

/*
 * DISCLAIMER:
 * The model contains (dimmz +2*HALO) * (dimmx +2HALO) * (dimmy +2*HALO) 
 * cells. For now, it assumes that there are no boundary conditions (sponge)
 * implemented.
 * The initial velocity model should be loaded taking into account this
 * criteria.
 */

int main(int argc, const char *argv[])
{
    if (argc != 2) {
        printf("Invalid argument! Schedule file is required.\n \
                Usage: %s <schedule_file>\n", argv[0]);
        abort();
    }

    if (argv[1] == NULL) {
        printf("Invalid argument! Schedule file path is NULL.\n \
                Usage: %s <schedule_file>\n", argv[0]);
        abort();
    }

    /* set seed for random number generator */
    srand(314);

    /* Load schedule file */
    schedule_t s = load_schedule(argv[1]);

    char foldername[500];
    char* fwipath = read_env_variable("FWIDIR");
    sprintf( foldername, "%s/data/inputmodels", fwipath);
    create_folder(foldername);


    /* Generate one velocity model per frequency */
    for(int i=0; i<s.nfreqs; i++)
    {
        real waveletFreq = s.freq[i];
        integer dimmz    = s.dimmz[i];
        integer dimmy    = s.dimmy[i];
        integer dimmx    = s.dimmx[i];

        print_info("Creating synthetic velocity input model for %f Hz freq", waveletFreq );

        /* generate complete path for output model */
        char modelname[500];
        sprintf( modelname, "%s/velocitymodel_%.2f.bin", foldername, waveletFreq );

        FILE* model = safe_fopen( modelname, "wb", __FILE__, __LINE__);

        /* compute number of cells per array */
        const integer cellsInVolume = dimmz * dimmy * dimmx;
        print_info("Number of cells in volume:"I, cellsInVolume);
        real *buffer = __malloc( ALIGN_REAL, sizeof(real) * cellsInVolume);

        /* safe dummy buffer */
        for(int i = 0; i < WRITTEN_FIELDS; i++)
        {
            set_array_to_random_real( buffer, cellsInVolume );
            safe_fwrite( buffer, sizeof(real), cellsInVolume, model, __FILE__, __LINE__);
        }

        /* free buffer */
        __free( buffer );

        /*  close model file */
        safe_fclose( modelname, model, __FILE__, __LINE__);

        print_info("Model %s created correctly", modelname);
    }

    schedule_free(s);

    print_info("End of the program");

    return 0;
}

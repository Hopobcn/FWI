/*
 * =====================================================================================
 *
 *       Filename:  GenerateInputModel.c
 *
 *    Description:  Generates input velocity model for the FWI code
 *
 *        Version:  1.0
 *        Created:  26/01/16 11:55:02
 *       Revision:  none
 *       Compiler:  icc
 *
 *         Author:  YOUR NAME (),
 *   Organization:
 *
 * =====================================================================================
 */

#include "fwi_kernel.h"


int main(int argc, const char *argv[])
{
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
        sprintf( modelname, "../InputModels/velocitymodel_%.2f.bin", waveletFreq );

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

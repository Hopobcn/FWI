/*
 * =====================================================================================
 *
 *       Filename:  fwi_kernel.c
 *
 *    Description:  kernel propagator implementation
 *
 *        Version:  1.0
 *        Created:  14/12/15 12:10:05
 *       Revision:  none
 *       Compiler:  icc
 *
 *         Author:  YOUR NAME (),
 *   Organization:
 *
 * =====================================================================================
 */

#include "fwi_kernel.h"

/*
 * Initializes and array of length "length" to a random number.
 */
void set_array_to_random_real( real* restrict array, const integer length)
{
    const real randvalue = rand() / (1.0 * RAND_MAX);

#ifdef DEBUG
    fprintf(stderr, "Array is being initialized to %f\n", randvalue);
#endif

    for( integer i = 0; i < length; i++ )
        array[i] = randvalue;
}

void check_memory_shot( const integer numberOfCells,
                        coeff_t *c,
                        s_t     *s,
                        v_t     *v,
                        real    *rho)
{
#ifdef DEBUG
    fprintf(stderr, "Checking memory shot values\n");

    real value;
    for( int i=0; i < numberOfCells; i++)
    {
        value = c->c11[i];
        value = c->c12[i];
        value = c->c13[i];
        value = c->c14[i];
        value = c->c15[i];
        value = c->c16[i];

        value = c->c22[i];
        value = c->c23[i];
        value = c->c24[i];
        value = c->c25[i];
        value = c->c26[i];

        value = c->c33[i];
        value = c->c34[i];
        value = c->c35[i];
        value = c->c36[i];

        value = c->c44[i];
        value = c->c45[i];
        value = c->c46[i];
        value =
        value = c->c55[i];
        value = c->c56[i];
        value = c->c66[i];

        value = v->tl.u[i];
        value = v->tl.v[i];
        value = v->tl.w[i];

        value = v->tr.u[i];
        value = v->tr.v[i];
        value = v->tr.w[i];

        value = v->bl.u[i];
        value = v->bl.v[i];
        value = v->bl.w[i];
        value =
        value = v->br.u[i];
        value = v->br.v[i];
        value = v->br.w[i];

        value = s->tl.zz[i];
        value = s->tl.xz[i];
        value = s->tl.yz[i];
        value = s->tl.xx[i];
        value = s->tl.xy[i];
        value = s->tl.yy[i];

        value = s->tr.zz[i];
        value = s->tr.xz[i];
        value = s->tr.yz[i];
        value = s->tr.xx[i];
        value = s->tr.xy[i];
        value = s->tr.yy[i];

        value = s->bl.zz[i];
        value = s->bl.xz[i];
        value = s->bl.yz[i];
        value = s->bl.xx[i];
        value = s->bl.xy[i];
        value = s->bl.yy[i];

        value = s->br.zz[i];
        value = s->br.xz[i];
        value = s->br.yz[i];
        value = s->br.xx[i];
        value = s->br.xy[i];
        value = s->br.yy[i];

        value = rho[i];
    }
#endif
};

void alloc_memory_shot( const integer numberOfCells,
                        coeff_t *c,
                        s_t     *s,
                        v_t     *v,
                        real    **rho)
{
    const integer size = numberOfCells * sizeof(real);

    fprintf(stderr, "ptr size = " I " bytes ("I" elements)\n", size, numberOfCells);

    /* allocate coefficients */
    c->c11 = (real*) __malloc( ALIGN_REAL, size);
    c->c12 = (real*) __malloc( ALIGN_REAL, size);
    c->c13 = (real*) __malloc( ALIGN_REAL, size);
    c->c14 = (real*) __malloc( ALIGN_REAL, size);
    c->c15 = (real*) __malloc( ALIGN_REAL, size);
    c->c16 = (real*) __malloc( ALIGN_REAL, size);

    c->c22 = (real*) __malloc( ALIGN_REAL, size);
    c->c23 = (real*) __malloc( ALIGN_REAL, size);
    c->c24 = (real*) __malloc( ALIGN_REAL, size);
    c->c25 = (real*) __malloc( ALIGN_REAL, size);
    c->c26 = (real*) __malloc( ALIGN_REAL, size);

    c->c33 = (real*) __malloc( ALIGN_REAL, size);
    c->c34 = (real*) __malloc( ALIGN_REAL, size);
    c->c35 = (real*) __malloc( ALIGN_REAL, size);
    c->c36 = (real*) __malloc( ALIGN_REAL, size);

    c->c44 = (real*) __malloc( ALIGN_REAL, size);
    c->c45 = (real*) __malloc( ALIGN_REAL, size);
    c->c46 = (real*) __malloc( ALIGN_REAL, size);

    c->c55 = (real*) __malloc( ALIGN_REAL, size);
    c->c56 = (real*) __malloc( ALIGN_REAL, size);
    c->c66 = (real*) __malloc( ALIGN_REAL, size);

    /* allocate velocity components */
    v->tl.u = (real*) __malloc( ALIGN_REAL, size);
    v->tl.v = (real*) __malloc( ALIGN_REAL, size);
    v->tl.w = (real*) __malloc( ALIGN_REAL, size);

    v->tr.u = (real*) __malloc( ALIGN_REAL, size);
    v->tr.v = (real*) __malloc( ALIGN_REAL, size);
    v->tr.w = (real*) __malloc( ALIGN_REAL, size);

    v->bl.u = (real*) __malloc( ALIGN_REAL, size);
    v->bl.v = (real*) __malloc( ALIGN_REAL, size);
    v->bl.w = (real*) __malloc( ALIGN_REAL, size);

    v->br.u = (real*) __malloc( ALIGN_REAL, size);
    v->br.v = (real*) __malloc( ALIGN_REAL, size);
    v->br.w = (real*) __malloc( ALIGN_REAL, size);

    /* allocate stress components   */
    s->tl.zz = (real*) __malloc( ALIGN_REAL, size);
    s->tl.xz = (real*) __malloc( ALIGN_REAL, size);
    s->tl.yz = (real*) __malloc( ALIGN_REAL, size);
    s->tl.xx = (real*) __malloc( ALIGN_REAL, size);
    s->tl.xy = (real*) __malloc( ALIGN_REAL, size);
    s->tl.yy = (real*) __malloc( ALIGN_REAL, size);

    s->tr.zz = (real*) __malloc( ALIGN_REAL, size);
    s->tr.xz = (real*) __malloc( ALIGN_REAL, size);
    s->tr.yz = (real*) __malloc( ALIGN_REAL, size);
    s->tr.xx = (real*) __malloc( ALIGN_REAL, size);
    s->tr.xy = (real*) __malloc( ALIGN_REAL, size);
    s->tr.yy = (real*) __malloc( ALIGN_REAL, size);

    s->bl.zz = (real*) __malloc( ALIGN_REAL, size);
    s->bl.xz = (real*) __malloc( ALIGN_REAL, size);
    s->bl.yz = (real*) __malloc( ALIGN_REAL, size);
    s->bl.xx = (real*) __malloc( ALIGN_REAL, size);
    s->bl.xy = (real*) __malloc( ALIGN_REAL, size);
    s->bl.yy = (real*) __malloc( ALIGN_REAL, size);

    s->br.zz = (real*) __malloc( ALIGN_REAL, size);
    s->br.xz = (real*) __malloc( ALIGN_REAL, size);
    s->br.yz = (real*) __malloc( ALIGN_REAL, size);
    s->br.xx = (real*) __malloc( ALIGN_REAL, size);
    s->br.xy = (real*) __malloc( ALIGN_REAL, size);
    s->br.yy = (real*) __malloc( ALIGN_REAL, size);

    /* allocate density array       */
    *rho = (real*) __malloc( ALIGN_REAL, size);
};

void free_memory_shot( coeff_t *c,
                       s_t     *s,
                       v_t     *v,
                       real    **rho)
{
    

    const int ngpus = acc_get_num_devices( acc_device_nvidia );
    
    #pragma omp parallel for schedule(static, 1)
    for (int gpu = 0; gpu < ngpus; gpu++)
    {
        

        acc_set_device_num(gpu, acc_device_nvidia);
        
        #pragma acc wait
        #pragma acc exit data delete(v->tl.u, v->tl.v, v->tl.w) \
                              delete(v->tr.u, v->tr.v, v->tr.w) \
                              delete(v->bl.u, v->bl.v, v->bl.w) \
                              delete(v->br.u, v->br.v, v->br.w) \
                              delete(s->tl.zz, s->tl.xz, s->tl.yz, s->tl.xx, s->tl.xy, s->tl.yy) \
                              delete(s->tr.zz, s->tr.xz, s->tr.yz, s->tr.xx, s->tr.xy, s->tr.yy) \
                              delete(s->bl.zz, s->bl.xz, s->bl.yz, s->bl.xx, s->bl.xy, s->bl.yy) \
                              delete(s->br.zz, s->br.xz, s->br.yz, s->br.xx, s->br.xy, s->br.yy) \
                              delete(c->c11, c->c12, c->c13, c->c14, c->c15, c->c16) \
                              delete(c->c22, c->c23, c->c24, c->c25, c->c26) \
                              delete(c->c33, c->c34, c->c35, c->c36) \
                              delete(c->c44, c->c45, c->c46) \
                              delete(c->c55, c->c56) \
                              delete(c->c66) \
                              delete(rho)
        
    }


    /* deallocate coefficients */
    __free( (void*) c->c11 );
    __free( (void*) c->c12 );
    __free( (void*) c->c13 );
    __free( (void*) c->c14 );
    __free( (void*) c->c15 );
    __free( (void*) c->c16 );

    __free( (void*) c->c22 );
    __free( (void*) c->c23 );
    __free( (void*) c->c24 );
    __free( (void*) c->c25 );
    __free( (void*) c->c26 );
    __free( (void*) c->c33 );

    __free( (void*) c->c34 );
    __free( (void*) c->c35 );
    __free( (void*) c->c36 );

    __free( (void*) c->c44 );
    __free( (void*) c->c45 );
    __free( (void*) c->c46 );

    __free( (void*) c->c55 );
    __free( (void*) c->c56 );

    __free( (void*) c->c66 );

    /* deallocate velocity components */
    __free( (void*) v->tl.u );
    __free( (void*) v->tl.v );
    __free( (void*) v->tl.w );

    __free( (void*) v->tr.u );
    __free( (void*) v->tr.v );
    __free( (void*) v->tr.w );

    __free( (void*) v->bl.u );
    __free( (void*) v->bl.v );
    __free( (void*) v->bl.w );

    __free( (void*) v->br.u );
    __free( (void*) v->br.v );
    __free( (void*) v->br.w );

    /* deallocate stres components   */
    __free( (void*) s->tl.zz );
    __free( (void*) s->tl.xz );
    __free( (void*) s->tl.yz );
    __free( (void*) s->tl.xx );
    __free( (void*) s->tl.xy );
    __free( (void*) s->tl.yy );

    __free( (void*) s->tr.zz );
    __free( (void*) s->tr.xz );
    __free( (void*) s->tr.yz );
    __free( (void*) s->tr.xx );
    __free( (void*) s->tr.xy );
    __free( (void*) s->tr.yy );

    __free( (void*) s->bl.zz );
    __free( (void*) s->bl.xz );
    __free( (void*) s->bl.yz );
    __free( (void*) s->bl.xx );
    __free( (void*) s->bl.xy );
    __free( (void*) s->bl.yy );

    __free( (void*) s->br.zz );
    __free( (void*) s->br.xz );
    __free( (void*) s->br.yz );
    __free( (void*) s->br.xx );
    __free( (void*) s->br.xy );
    __free( (void*) s->br.yy );


    /* deallocate density array       */
    __free( (void*) *rho );

    
};

/*
 * Loads initial values from coeffs, stress and velocity.
 */
 void load_initial_model ( const real    waveletFreq,
                           const integer numberOfCells,
                           coeff_t *c,
                           s_t     *s,
                           v_t     *v,
                           real    *rho)
{
    const integer size = numberOfCells * sizeof(real);
    
       /* initialize coefficients */
    set_array_to_random_real( c->c11, numberOfCells);
    set_array_to_random_real( c->c12, numberOfCells);
    set_array_to_random_real( c->c13, numberOfCells);
    set_array_to_random_real( c->c14, numberOfCells);
    set_array_to_random_real( c->c15, numberOfCells);
    set_array_to_random_real( c->c16, numberOfCells);
    set_array_to_random_real( c->c22, numberOfCells);
    set_array_to_random_real( c->c23, numberOfCells);
    set_array_to_random_real( c->c24, numberOfCells);
    set_array_to_random_real( c->c25, numberOfCells);
    set_array_to_random_real( c->c26, numberOfCells);
    set_array_to_random_real( c->c33, numberOfCells);
    set_array_to_random_real( c->c34, numberOfCells);
    set_array_to_random_real( c->c35, numberOfCells);
    set_array_to_random_real( c->c36, numberOfCells);
    set_array_to_random_real( c->c44, numberOfCells);
    set_array_to_random_real( c->c45, numberOfCells);
    set_array_to_random_real( c->c46, numberOfCells);
    set_array_to_random_real( c->c55, numberOfCells);
    set_array_to_random_real( c->c56, numberOfCells);
    set_array_to_random_real( c->c66, numberOfCells);
    
    /* initialize stress */
    memset( s->tl.zz, 0, size);
    memset( s->tl.xz, 0, size);
    memset( s->tl.yz, 0, size);
    memset( s->tl.xx, 0, size);
    memset( s->tl.xy, 0, size);
    memset( s->tl.yy, 0, size);
    memset( s->tr.zz, 0, size);
    memset( s->tr.xz, 0, size);
    memset( s->tr.yz, 0, size);
    memset( s->tr.xx, 0, size);
    memset( s->tr.xy, 0, size);
    memset( s->tr.yy, 0, size);
    memset( s->bl.zz, 0, size);
    memset( s->bl.xz, 0, size);
    memset( s->bl.yz, 0, size);
    memset( s->bl.xx, 0, size);
    memset( s->bl.xy, 0, size);
    memset( s->bl.yy, 0, size);
    memset( s->br.zz, 0, size);
    memset( s->br.xz, 0, size);
    memset( s->br.yz, 0, size);
    memset( s->br.xx, 0, size);
    memset( s->br.xy, 0, size);
    memset( s->br.yy, 0, size);

#ifdef DO_NOT_PERFORM_IO /* initalize velocity components */

    set_array_to_random_real( v->tl.u, numberOfCells );
    set_array_to_random_real( v->tl.v, numberOfCells );
    set_array_to_random_real( v->tl.w, numberOfCells );
    set_array_to_random_real( v->tr.u, numberOfCells );
    set_array_to_random_real( v->tr.v, numberOfCells );
    set_array_to_random_real( v->tr.w, numberOfCells );
    set_array_to_random_real( v->bl.u, numberOfCells );
    set_array_to_random_real( v->bl.v, numberOfCells );
    set_array_to_random_real( v->bl.w, numberOfCells );
    set_array_to_random_real( v->br.u, numberOfCells );
    set_array_to_random_real( v->br.v, numberOfCells );
    set_array_to_random_real( v->br.w, numberOfCells );
    
    /* initialize rho */
    set_array_to_random_real( rho, numberOfCells );

#else /* load velocity model from external file */

    /* open initial model, binary file */
    char modelname[300];
    sprintf( modelname, "../InputModels/velocitymodel_%.2f.bin", waveletFreq );

    fprintf(stderr, "Loading input model %s from disk (this could take a while)\n", modelname);

    FILE* model = safe_fopen( modelname, "rb", __FILE__, __LINE__ );

    /* initalize velocity components */
    safe_fread( v->tl.u, sizeof(real), numberOfCells, model, __FILE__, __LINE__ );
    safe_fread( v->tl.v, sizeof(real), numberOfCells, model, __FILE__, __LINE__ );
    safe_fread( v->tl.w, sizeof(real), numberOfCells, model, __FILE__, __LINE__ );
    safe_fread( v->tr.u, sizeof(real), numberOfCells, model, __FILE__, __LINE__ );
    safe_fread( v->tr.v, sizeof(real), numberOfCells, model, __FILE__, __LINE__ );
    safe_fread( v->tr.w, sizeof(real), numberOfCells, model, __FILE__, __LINE__ );
    safe_fread( v->bl.u, sizeof(real), numberOfCells, model, __FILE__, __LINE__ );
    safe_fread( v->bl.v, sizeof(real), numberOfCells, model, __FILE__, __LINE__ );
    safe_fread( v->bl.w, sizeof(real), numberOfCells, model, __FILE__, __LINE__ );
    safe_fread( v->br.u, sizeof(real), numberOfCells, model, __FILE__, __LINE__ );
    safe_fread( v->br.v, sizeof(real), numberOfCells, model, __FILE__, __LINE__ );
    safe_fread( v->br.w, sizeof(real), numberOfCells, model, __FILE__, __LINE__ );

    /* close model file */
    safe_fclose ( "velocitymodel.bin", model, __FILE__, __LINE__ );

#endif /* end of DDO_NOT_PERFORM_IO clause */

    const integer datalen = numberOfCells;

    const real* cc11 = c->c11;
    const real* cc12 = c->c12;
    const real* cc13 = c->c13;
    const real* cc14 = c->c14;
    const real* cc15 = c->c15;
    const real* cc16 = c->c16;
                             
    const real* cc22 = c->c22;
    const real* cc23 = c->c23;
    const real* cc24 = c->c24;
    const real* cc25 = c->c25;
    const real* cc26 = c->c26;
                             
    const real* cc33 = c->c33;
    const real* cc34 = c->c34;
    const real* cc35 = c->c35;
    const real* cc36 = c->c36;
                             
    const real* cc44 = c->c44;
    const real* cc45 = c->c45;
    const real* cc46 = c->c46;
                             
    const real* cc55 = c->c55;
    const real* cc56 = c->c56;
                       
    const real* cc66 = c->c66;
          
    const real* vtlu = v->tl.u;
    const real* vtlv = v->tl.v;
    const real* vtlw = v->tl.w;
                              
    const real* vtru = v->tr.u;
    const real* vtrv = v->tr.v;
    const real* vtrw = v->tr.w;
                              
    const real* vblu = v->bl.u;
    const real* vblv = v->bl.v;
    const real* vblw = v->bl.w;
                              
    const real* vbru = v->br.u;
    const real* vbrv = v->br.v;
    const real* vbrw = v->br.w;
       
    const real* stlzz = s->tl.zz;
    const real* stlxz = s->tl.xz;
    const real* stlyz = s->tl.yz;
    const real* stlxx = s->tl.xx;
    const real* stlxy = s->tl.xy;
    const real* stlyy = s->tl.yy;
                                
    const real* strzz = s->tr.zz;
    const real* strxz = s->tr.xz;
    const real* stryz = s->tr.yz;
    const real* strxx = s->tr.xx;
    const real* strxy = s->tr.xy;
    const real* stryy = s->tr.yy;
                                
    const real* sblzz = s->bl.zz;
    const real* sblxz = s->bl.xz;
    const real* sblyz = s->bl.yz;
    const real* sblxx = s->bl.xx;
    const real* sblxy = s->bl.xy;
    const real* sblyy = s->bl.yy;
                                
    const real* sbrzz = s->br.zz;
    const real* sbrxz = s->br.xz;
    const real* sbryz = s->br.yz;
    const real* sbrxx = s->br.xx;
    const real* sbrxy = s->br.xy;
    const real* sbryy = s->br.yy;

    const real* rrho  = rho;

    const int ngpus = acc_get_num_devices( acc_device_nvidia );

    //#pragma omp parallel for schedule(static, 1)    
    for (int gpu = 0; gpu < ngpus; gpu++)
    {
        

        acc_set_device_num(gpu, acc_device_nvidia);
       
        #pragma acc enter data copyin(vtlu[0:datalen], vtlv[0:datalen], vtlw[0:datalen]) \
                               copyin(vtru[0:datalen], vtrv[0:datalen], vtrw[0:datalen]) \
                               copyin(vblu[0:datalen], vblv[0:datalen], vblw[0:datalen]) \
                               copyin(vbru[0:datalen], vbrv[0:datalen], vbrw[0:datalen]) \
                               create(stlzz[0:datalen], stlxz[0:datalen], stlyz[0:datalen], stlxx[0:datalen], stlxy[0:datalen], stlyy[0:datalen]) \
                               create(strzz[0:datalen], strxz[0:datalen], stryz[0:datalen], strxx[0:datalen], strxy[0:datalen], stryy[0:datalen]) \
                               create(sblzz[0:datalen], sblxz[0:datalen], sblyz[0:datalen], sblxx[0:datalen], sblxy[0:datalen], sblyy[0:datalen]) \
                               create(sbrzz[0:datalen], sbrxz[0:datalen], sbryz[0:datalen], sbrxx[0:datalen], sbrxy[0:datalen], sbryy[0:datalen]) \
                               copyin(cc11[0:datalen], cc12[0:datalen], cc13[0:datalen], cc14[0:datalen], cc15[0:datalen], cc16[0:datalen]) \
                               copyin(cc22[0:datalen], cc23[0:datalen], cc24[0:datalen], cc25[0:datalen], cc26[0:datalen]) \
                               copyin(cc33[0:datalen], cc34[0:datalen], cc35[0:datalen], cc36[0:datalen]) \
                               copyin(cc44[0:datalen], cc45[0:datalen], cc46[0:datalen]) \
                               copyin(cc55[0:datalen], cc56[0:datalen]) \
                               copyin(cc66[0:datalen]) \
                               copyin(rrho[0:datalen]) 
       
        #pragma acc wait // wait for the copies to be finished before start executing kernels (we could optimize this)
        
        
    }

};


/*
 * Saves the complete velocity field to disk.
 */
void write_snapshot(char *folder,
                    int suffix,
                    v_t *v,
                    const integer numberOfCells)
{
#ifdef DO_NOT_PERFORM_IO
    fprintf(stderr, "Warning: We are not doing any IO here (%s)\n", __FUNCTION__);

#else
    #pragma acc update self(v->tr.u[0:numberOfCells], v->tr.v[0:numberOfCells], v->tr.w[0:numberOfCells]) \
                       self(v->tl.u[0:numberOfCells], v->tl.v[0:numberOfCells], v->tl.w[0:numberOfCells]) \
                       self(v->br.u[0:numberOfCells], v->br.v[0:numberOfCells], v->br.w[0:numberOfCells]) \
                       self(v->bl.u[0:numberOfCells], v->bl.v[0:numberOfCells], v->bl.w[0:numberOfCells])


    char fname[300];
    sprintf(fname,"%s/snapshot.%05d.bin", folder, suffix);

    FILE *snapshot = safe_fopen(fname,"wb", __FILE__, __LINE__ );

    safe_fwrite( v->tr.u, sizeof(real), numberOfCells, snapshot, __FILE__, __LINE__ );
    safe_fwrite( v->tr.v, sizeof(real), numberOfCells, snapshot, __FILE__, __LINE__ );
    safe_fwrite( v->tr.w, sizeof(real), numberOfCells, snapshot, __FILE__, __LINE__ );
    
    safe_fwrite( v->tl.u, sizeof(real), numberOfCells, snapshot, __FILE__, __LINE__ );
    safe_fwrite( v->tl.v, sizeof(real), numberOfCells, snapshot, __FILE__, __LINE__ );
    safe_fwrite( v->tl.w, sizeof(real), numberOfCells, snapshot, __FILE__, __LINE__ );

    safe_fwrite( v->br.u, sizeof(real), numberOfCells, snapshot, __FILE__, __LINE__ );
    safe_fwrite( v->br.v, sizeof(real), numberOfCells, snapshot, __FILE__, __LINE__ );
    safe_fwrite( v->br.w, sizeof(real), numberOfCells, snapshot, __FILE__, __LINE__ );
    
    safe_fwrite( v->bl.u, sizeof(real), numberOfCells, snapshot, __FILE__, __LINE__ );
    safe_fwrite( v->bl.v, sizeof(real), numberOfCells, snapshot, __FILE__, __LINE__ );
    safe_fwrite( v->bl.w, sizeof(real), numberOfCells, snapshot, __FILE__, __LINE__ );

    if ( fclose(snapshot)!=0)
        fprintf(stderr,"Error closing file %s\n", fname);

#endif
};

/*
 * Reads the complete velocity field from disk.
 */
void read_snapshot(char *folder,
                   int suffix,
                   v_t *v,
                   const integer numberOfCells)
{
#ifdef DO_NOT_PERFORM_IO
    fprintf(stderr, "Warning: We are not doing any IO here (%s)\n", __FUNCTION__);

#else
    char fname[300];
    sprintf(fname,"%s/snapshot.%05d.bin", folder, suffix);

    FILE *snapshot = safe_fopen(fname,"rb", __FILE__, __LINE__ );

    safe_fread( v->tr.u, sizeof(real), numberOfCells, snapshot, __FILE__, __LINE__ );
    safe_fread( v->tr.v, sizeof(real), numberOfCells, snapshot, __FILE__, __LINE__ );
    safe_fread( v->tr.w, sizeof(real), numberOfCells, snapshot, __FILE__, __LINE__ );
    
    safe_fread( v->tl.u, sizeof(real), numberOfCells, snapshot, __FILE__, __LINE__ );
    safe_fread( v->tl.v, sizeof(real), numberOfCells, snapshot, __FILE__, __LINE__ );
    safe_fread( v->tl.w, sizeof(real), numberOfCells, snapshot, __FILE__, __LINE__ );

    safe_fread( v->br.u, sizeof(real), numberOfCells, snapshot, __FILE__, __LINE__ );
    safe_fread( v->br.v, sizeof(real), numberOfCells, snapshot, __FILE__, __LINE__ );
    safe_fread( v->br.w, sizeof(real), numberOfCells, snapshot, __FILE__, __LINE__ );
    
    safe_fread( v->bl.u, sizeof(real), numberOfCells, snapshot, __FILE__, __LINE__ );
    safe_fread( v->bl.v, sizeof(real), numberOfCells, snapshot, __FILE__, __LINE__ );
    safe_fread( v->bl.w, sizeof(real), numberOfCells, snapshot, __FILE__, __LINE__ );

    if ( fclose(snapshot)!=0 )
        fprintf(stderr,"Error closing file %s\n", fname);

    #pragma acc update device(v->tr.u[0:numberOfCells], v->tr.v[0:numberOfCells], v->tr.w[0:numberOfCells]) \
                       device(v->tl.u[0:numberOfCells], v->tl.v[0:numberOfCells], v->tl.w[0:numberOfCells]) \
                       device(v->br.u[0:numberOfCells], v->br.v[0:numberOfCells], v->br.w[0:numberOfCells]) \
                       device(v->bl.u[0:numberOfCells], v->bl.v[0:numberOfCells], v->bl.w[0:numberOfCells])

#endif
};

void propagate_shot (time_d        direction,
                     v_t           v,
                     s_t           s,
                     coeff_t       coeffs,
                     real          *rho,
                     int           timesteps,
                     int           ntbwd,
                     real          dt,
                     real          dzi,
                     real          dxi,
                     real          dyi,
                     integer       nz0,
                     integer       nzf,
                     integer       nx0,
                     integer       nxf,
                     integer       ny0,
                     integer       nyf,
                     integer       stacki,
                     char          *folder,
                     real          *dataflush,
                     integer       datalen,
                     integer       dimmz,
                     integer       dimmx)
{
    for(int t=0; t < timesteps; t++)
    {
        fprintf(stderr, "Computing %d-th timestep\n", t);

        /* perform IO */
        if ( t%stacki == 0 && direction == BACKWARD) read_snapshot(folder, ntbwd-t, &v, datalen);

        const int ngpus              = acc_get_num_devices( acc_device_nvidia );
        const int num_planes_per_gpu = (nyf-2*HALO) / ngpus;
        const int remainder          = (nyf-2*HALO) % ngpus;
        
        #pragma omp parallel for schedule(static, 1)
        for (int gpu = 0; gpu < ngpus; gpu++)
        {
            const int roffset = ((gpu < remainder) ? 1 : 0);
            const int accumof = ( remainder ); // TODO: control situations when remainder>0!
            const int ny0i    = ny0  + gpu * (num_planes_per_gpu + accumof) + ((gpu > 0) ? HALO : 0);
            const int nyfi    = ny0i +       (num_planes_per_gpu + roffset) + (((gpu > 0) && (gpu < ngpus-1)) ? 0 : HALO);
            const int plane_size = ((nzf-nz0) + 2*HALO) * ((nxf-nx0) + 2*HALO);
            
            //fprintf(stderr, "gpu %d ny0i %d nybi %d -- (ny0 %d nyf %d) omp_thread_id %d\n", 
            //        gpu, ny0i, nyfi, ny0, nyf, omp_get_thread_num());

            acc_set_device_num(gpu, acc_device_nvidia);
 
            /* ------------------------------------------------------------------------------ */
            /*                      VELOCITY COMPUTATION                                      */
            /* ------------------------------------------------------------------------------ */
           
            if (gpu > 0) {
                /* Phase 1. Computation of the left-most planes of the domain */
                velocity_propagator(v, s, coeffs, rho, dt, dzi, dxi, dyi,
                                    nz0 +   HALO,
                                    nzf -   HALO,
                                    nx0 +   HALO,
                                    nxf -   HALO,
                                    ny0i,
                                    ny0i+   HALO,
                                    dimmz, dimmx,
                                    ONE_L); 
            } 
            if (gpu < ngpus - 1) {
                /* Phase 1. Computation of the right-most planes of the domain */
                velocity_propagator(v, s, coeffs, rho, dt, dzi, dxi, dyi,
                                    nz0 +   HALO,
                                    nzf -   HALO,
                                    nx0 +   HALO,
                                    nxf -   HALO,
                                    nyfi-   HALO,
                                    nyfi,
                                    dimmz, dimmx,
                                    ONE_R);
            } 
            /* Phase 2. Computation of the central planes. */
            velocity_propagator(v, s, coeffs, rho, dt, dzi, dxi, dyi,
                                nz0 +   HALO,
                                nzf -   HALO,
                                nx0 +   HALO,
                                nxf -   HALO,
                                ny0i+   HALO,
                                nyfi-   HALO,
                                dimmz, dimmx,
                                TWO);

            /* Boundary exchange for velocity values */
            exchange_velocity_boundaries( v, plane_size, gpu, ngpus, nyfi, ny0i);

            #pragma acc wait(ONE_L, ONE_R, TWO)

            /* ------------------------------------------------------------------------------ */
            /*                        STRESS COMPUTATION                                      */
            /* ------------------------------------------------------------------------------ */

            if (gpu > 0) {
                /* Phase 1. Computation of the left-most planes of the domain */
                stress_propagator ( s, v, coeffs, rho, dt, dzi, dxi, dyi, 
                                    nz0 +   HALO,
                                    nzf -   HALO,
                                    nx0 +   HALO,
                                    nxf -   HALO,
                                    ny0i,
                                    ny0i+   HALO,
                                    dimmz, dimmx,
                                    ONE_L);
            }
            if (gpu < ngpus - 1) {
                /* Phase 1. Computation of the right-most planes of the domain */
                stress_propagator ( s, v, coeffs, rho, dt, dzi, dxi, dyi, 
                                    nz0 +   HALO,
                                    nzf -   HALO,
                                    nx0 +   HALO,
                                    nxf -   HALO,
                                    nyfi-   HALO,
                                    nyfi,
                                    dimmz, dimmx,
                                    ONE_R);
            }

            /* Phase 2 computation. Central planes of the domain */
            stress_propagator ( s, v, coeffs, rho, dt, dzi, dxi, dyi, 
                                nz0 +   HALO,
                                nzf -   HALO,
                                nx0 +   HALO,
                                nxf -   HALO,
                                ny0i+   HALO,
                                nyfi-   HALO,
                                dimmz, dimmx,
                                TWO);

            /* Boundary exchange for stress values */
            exchange_stress_boundaries( s, plane_size, gpu, ngpus, nyfi, ny0i);

            #pragma acc wait(ONE_L, ONE_R, TWO, H2D, D2H)
        } 
              /* perform IO */
        if ( t%stacki == 0 && direction == FORWARD) write_snapshot(folder, ntbwd-t, &v, datalen);
    }
};

/*
NAME:exchange_boundaries
PURPOSE: data exchanges between the boundary layers of the analyzed volume

v                   (in) struct containing velocity arrays (4 points / cell x 3 components / point = 12 arrays)
plane_size          (in) Number of elements per plane to be exchanged
gpu                 (in) gpu id
ngpus               (in) number of gpus
nyf                 (in) final plane to be exchanged
ny0                 (in) intial plane to be exchanged

RETURN none
*/
void exchange_velocity_boundaries ( v_t v, 
                                    const integer plane_size, 
                                    const integer gpu,
                                    const integer ngpus,
                                    const integer nyf, 
                                    const integer ny0 )
{
    const integer num_planes = (nyf-ny0);
    const integer nelems     = num_planes * plane_size;

    if ( gpu != 0 )
    {
        acc_set_device_num(gpu-1, acc_device_nvidia);
 
        #pragma acc update device(v.tl.u[ny0:nelems]) wait(D2H) asynch(H2D)
        #pragma acc update device(v.tl.v[ny0:nelems]) wait(D2H) asynch(H2D)
        #pragma acc update device(v.tl.w[ny0:nelems]) wait(D2H) asynch(H2D)
 
        #pragma acc update device(v.tr.u[ny0:nelems]) wait(D2H) asynch(H2D)
        #pragma acc update device(v.tr.v[ny0:nelems]) wait(D2H) asynch(H2D)
        #pragma acc update device(v.tr.w[ny0:nelems]) wait(D2H) asynch(H2D)
       
        #pragma acc update device(v.bl.u[ny0:nelems]) wait(D2H) asynch(H2D)
        #pragma acc update device(v.bl.v[ny0:nelems]) wait(D2H) asynch(H2D)
        #pragma acc update device(v.bl.w[ny0:nelems]) wait(D2H) asynch(H2D)
       
        #pragma acc update device(v.br.u[ny0:nelems]) wait(D2H) asynch(H2D)
        #pragma acc update device(v.br.v[ny0:nelems]) wait(D2H) asynch(H2D)
        #pragma acc update device(v.br.w[ny0:nelems]) wait(D2H) asynch(H2D)
    }

    if ( gpu != ngpus-1 )
    {
        acc_set_device_num(gpu+1, acc_device_nvidia);
 
        #pragma acc update device(v.tl.u[nyf:nelems]) wait(D2H) asynch(H2D)
        #pragma acc update device(v.tl.v[nyf:nelems]) wait(D2H) asynch(H2D)
        #pragma acc update device(v.tl.w[nyf:nelems]) wait(D2H) asynch(H2D)
 
        #pragma acc update device(v.tr.u[nyf:nelems]) wait(D2H) asynch(H2D)
        #pragma acc update device(v.tr.v[nyf:nelems]) wait(D2H) asynch(H2D)
        #pragma acc update device(v.tr.w[nyf:nelems]) wait(D2H) asynch(H2D)
       
        #pragma acc update device(v.bl.u[nyf:nelems]) wait(D2H) asynch(H2D)
        #pragma acc update device(v.bl.v[nyf:nelems]) wait(D2H) asynch(H2D)
        #pragma acc update device(v.bl.w[nyf:nelems]) wait(D2H) asynch(H2D)
       
        #pragma acc update device(v.br.u[nyf:nelems]) wait(D2H) asynch(H2D)
        #pragma acc update device(v.br.v[nyf:nelems]) wait(D2H) asynch(H2D)
        #pragma acc update device(v.br.w[nyf:nelems]) wait(D2H) asynch(H2D)
    }
    acc_set_device_num(gpu, acc_device_nvidia);
};

/*
NAME:exchange_stress_boundaries
PURPOSE: data exchanges between the boundary layers of the analyzed volume

s                   (in) struct containing stress arrays (4 points / cell x 6 components / point = 24 arrays)
numElement          (in) Number of elements to exchange
idxt                (in) identifier related to the folder
nyf                 (in) final plane to be exchanged
ny0                 (in) intial plane to be exchanged

RETURN none
*/
void exchange_stress_boundaries ( s_t s, 
                                  const integer plane_size, 
                                  const integer gpu,
                                  const integer ngpus,
                                  const integer nyf, 
                                  const integer ny0 )
{
    const integer nplanes = HALO;
    const integer nelems  = nplanes * plane_size;

    if ( gpu != 0 ) 
    {

        acc_set_device_num(gpu-1, acc_device_nvidia);
 
        #pragma acc update device(s.tl.zz[ny0:nelems]) wait(D2H) asynch(H2D)
        #pragma acc update device(s.tl.xz[ny0:nelems]) wait(D2H) asynch(H2D)
        #pragma acc update device(s.tl.yz[ny0:nelems]) wait(D2H) asynch(H2D)
        #pragma acc update device(s.tl.xx[ny0:nelems]) wait(D2H) asynch(H2D)
        #pragma acc update device(s.tl.xy[ny0:nelems]) wait(D2H) asynch(H2D)
        #pragma acc update device(s.tl.yy[ny0:nelems]) wait(D2H) asynch(H2D)

        #pragma acc update device(s.tr.zz[ny0:nelems]) wait(D2H) asynch(H2D)
        #pragma acc update device(s.tr.xz[ny0:nelems]) wait(D2H) asynch(H2D)
        #pragma acc update device(s.tr.yz[ny0:nelems]) wait(D2H) asynch(H2D)
        #pragma acc update device(s.tr.xx[ny0:nelems]) wait(D2H) asynch(H2D)
        #pragma acc update device(s.tr.xy[ny0:nelems]) wait(D2H) asynch(H2D)
        #pragma acc update device(s.tr.yy[ny0:nelems]) wait(D2H) asynch(H2D)

        #pragma acc update device(s.bl.zz[ny0:nelems]) wait(D2H) asynch(H2D)
        #pragma acc update device(s.bl.xz[ny0:nelems]) wait(D2H) asynch(H2D)
        #pragma acc update device(s.bl.yz[ny0:nelems]) wait(D2H) asynch(H2D)
        #pragma acc update device(s.bl.xx[ny0:nelems]) wait(D2H) asynch(H2D)
        #pragma acc update device(s.bl.xy[ny0:nelems]) wait(D2H) asynch(H2D)
        #pragma acc update device(s.bl.yy[ny0:nelems]) wait(D2H) asynch(H2D)

        #pragma acc update device(s.br.zz[ny0:nelems]) wait(D2H) asynch(H2D)
        #pragma acc update device(s.br.xz[ny0:nelems]) wait(D2H) asynch(H2D)
        #pragma acc update device(s.br.yz[ny0:nelems]) wait(D2H) asynch(H2D)
        #pragma acc update device(s.br.xx[ny0:nelems]) wait(D2H) asynch(H2D)
        #pragma acc update device(s.br.xy[ny0:nelems]) wait(D2H) asynch(H2D)
        #pragma acc update device(s.br.yy[ny0:nelems]) wait(D2H) asynch(H2D)
    }

    if ( gpu != ngpus -1 )  
    {
        #pragma acc update device(s.tl.zz[nyf:nelems]) wait(D2H) asynch(H2D)
        #pragma acc update device(s.tl.xz[nyf:nelems]) wait(D2H) asynch(H2D)
        #pragma acc update device(s.tl.yz[nyf:nelems]) wait(D2H) asynch(H2D)
        #pragma acc update device(s.tl.xx[nyf:nelems]) wait(D2H) asynch(H2D)
        #pragma acc update device(s.tl.xy[nyf:nelems]) wait(D2H) asynch(H2D)
        #pragma acc update device(s.tl.yy[nyf:nelems]) wait(D2H) asynch(H2D)

        #pragma acc update device(s.tr.zz[nyf:nelems]) wait(D2H) asynch(H2D)
        #pragma acc update device(s.tr.xz[nyf:nelems]) wait(D2H) asynch(H2D)
        #pragma acc update device(s.tr.yz[nyf:nelems]) wait(D2H) asynch(H2D)
        #pragma acc update device(s.tr.xx[nyf:nelems]) wait(D2H) asynch(H2D)
        #pragma acc update device(s.tr.xy[nyf:nelems]) wait(D2H) asynch(H2D)
        #pragma acc update device(s.tr.yy[nyf:nelems]) wait(D2H) asynch(H2D)

        #pragma acc update device(s.bl.zz[nyf:nelems]) wait(D2H) asynch(H2D)
        #pragma acc update device(s.bl.xz[nyf:nelems]) wait(D2H) asynch(H2D)
        #pragma acc update device(s.bl.yz[nyf:nelems]) wait(D2H) asynch(H2D)
        #pragma acc update device(s.bl.xx[nyf:nelems]) wait(D2H) asynch(H2D)
        #pragma acc update device(s.bl.xy[nyf:nelems]) wait(D2H) asynch(H2D)
        #pragma acc update device(s.bl.yy[nyf:nelems]) wait(D2H) asynch(H2D)

        #pragma acc update device(s.br.zz[nyf:nelems]) wait(D2H) asynch(H2D)
        #pragma acc update device(s.br.xz[nyf:nelems]) wait(D2H) asynch(H2D)
        #pragma acc update device(s.br.yz[nyf:nelems]) wait(D2H) asynch(H2D)
        #pragma acc update device(s.br.xx[nyf:nelems]) wait(D2H) asynch(H2D)
        #pragma acc update device(s.br.xy[nyf:nelems]) wait(D2H) asynch(H2D)
        #pragma acc update device(s.br.yy[nyf:nelems]) wait(D2H) asynch(H2D)
    }
    acc_set_device_num(gpu, acc_device_nvidia);
};


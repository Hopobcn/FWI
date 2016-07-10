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

void check_memory_shot( const integer numberOfCells,
                        coeff_t *c,
                        s_t     *s,
                        v_t     *v,
                        real    *rho)
{
    log_info ( "Checking memory shot values");
            
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
    
    log_info ( "Checking memory shot values: OK");
};

void alloc_memory_shot( const integer numberOfCells,
                        coeff_t *c,
                        s_t     *s,
                        v_t     *v,
                        real    **rho)
{
    const integer size = numberOfCells * sizeof(real);
    
    log_info ( "ptr size = " I " bytes ("I" elements)", size, numberOfCells);
    
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

    log_info ("Allocated memory for shot successfully");
};

void free_memory_shot( coeff_t *c,
                       s_t     *s,
                       v_t     *v,
                       real    **rho)
{
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

    log_info ("Shot memory freed");
};

/*
 * Loads initial values from coeffs, stress and velocity.
 */
void load_initial_model ( const real    waveletFreq,
                          const integer dimmz,
                          const integer dimmx,
                          const integer dimmy,
                          coeff_t       *c,
                          s_t           *s,
                          v_t           *v,
                          real          *rho)
{
    /* primero tengo que situarme en mi mundo MPI */
    int rankLocal, rankSize;
    MPI_Comm_rank( MPI_COMM_WORLD, &rankLocal );
    MPI_Comm_size( MPI_COMM_WORLD, &rankSize  );

    log_info("Loading initial velocity model");

    /* Number of cells for HALO planes and computational volume */
    const integer cellsInHalo   = (dimmz + 2 * HALO ) * (dimmx + 2*HALO);
    const integer cellsInVolume = (dimmz + 2 * HALO ) * (dimmx + 2*HALO) * ( dimmz/rankSize );  
    const integer cellsInArray  = (dimmz + 2 * HALO ) + (dimmx + 2*HALO) * ((dimmz/rankSize) + 2*HALO); 


    /* Size in bytes for HALO planes and computational volume   */
    const integer UNUSED(bytesForHalo)     = cellsInHalo   * sizeof(real);
    const integer UNUSED(bytesForVolume)   = cellsInVolume * sizeof(real);
    
    /* initialize coefficients */
    memset( c->c11, 0, cellsInArray );
    memset( c->c12, 0, cellsInArray );
    memset( c->c13, 0, cellsInArray );
    memset( c->c14, 0, cellsInArray );
    memset( c->c15, 0, cellsInArray );
    memset( c->c16, 0, cellsInArray );
    memset( c->c22, 0, cellsInArray );
    memset( c->c23, 0, cellsInArray );
    memset( c->c24, 0, cellsInArray );
    memset( c->c25, 0, cellsInArray );
    memset( c->c26, 0, cellsInArray );
    memset( c->c33, 0, cellsInArray );
    memset( c->c34, 0, cellsInArray );
    memset( c->c35, 0, cellsInArray );
    memset( c->c36, 0, cellsInArray );
    memset( c->c44, 0, cellsInArray );
    memset( c->c45, 0, cellsInArray );
    memset( c->c46, 0, cellsInArray );
    memset( c->c55, 0, cellsInArray );
    memset( c->c56, 0, cellsInArray );
    memset( c->c66, 0, cellsInArray );
      
    /* initialize stress */
    memset( s->tl.zz, 0, cellsInArray );
    memset( s->tl.xz, 0, cellsInArray );
    memset( s->tl.yz, 0, cellsInArray );
    memset( s->tl.xx, 0, cellsInArray );
    memset( s->tl.xy, 0, cellsInArray );
    memset( s->tl.yy, 0, cellsInArray );
    memset( s->tr.zz, 0, cellsInArray );
    memset( s->tr.xz, 0, cellsInArray );
    memset( s->tr.yz, 0, cellsInArray );
    memset( s->tr.xx, 0, cellsInArray );
    memset( s->tr.xy, 0, cellsInArray );
    memset( s->tr.yy, 0, cellsInArray );
    memset( s->bl.zz, 0, cellsInArray );
    memset( s->bl.xz, 0, cellsInArray );
    memset( s->bl.yz, 0, cellsInArray );
    memset( s->bl.xx, 0, cellsInArray );
    memset( s->bl.xy, 0, cellsInArray );
    memset( s->bl.yy, 0, cellsInArray );
    memset( s->br.zz, 0, cellsInArray );
    memset( s->br.xz, 0, cellsInArray );
    memset( s->br.yz, 0, cellsInArray );
    memset( s->br.xx, 0, cellsInArray );
    memset( s->br.xy, 0, cellsInArray );
    memset( s->br.yy, 0, cellsInArray );

#ifndef DO_NOT_PERFOM_IO 

    /* open initial model, binary file */
    char modelname[300];    
    sprintf( modelname, "../InputModels/velocitymodel_%.2f.bin", waveletFreq );

    log_info ("Loading input model %s from disk (this could take a while)", modelname);

    FILE* model = safe_fopen( modelname, "rb", __FILE__, __LINE__ );

    /* seek to the correct position given by HALO planes */
    fseek ( model, bytesForVolume * rankLocal, SEEK_SET);
    
    /* initalize velocity components */
    safe_fread( v->tl.u, sizeof(real), bytesForVolume, model, __FILE__, __LINE__ );
    safe_fread( v->tl.v, sizeof(real), bytesForVolume, model, __FILE__, __LINE__ );
    safe_fread( v->tl.w, sizeof(real), bytesForVolume, model, __FILE__, __LINE__ );
    safe_fread( v->tr.u, sizeof(real), bytesForVolume, model, __FILE__, __LINE__ );
    safe_fread( v->tr.v, sizeof(real), bytesForVolume, model, __FILE__, __LINE__ );
    safe_fread( v->tr.w, sizeof(real), bytesForVolume, model, __FILE__, __LINE__ );
    safe_fread( v->bl.u, sizeof(real), bytesForVolume, model, __FILE__, __LINE__ );
    safe_fread( v->bl.v, sizeof(real), bytesForVolume, model, __FILE__, __LINE__ );
    safe_fread( v->bl.w, sizeof(real), bytesForVolume, model, __FILE__, __LINE__ );
    safe_fread( v->br.u, sizeof(real), bytesForVolume, model, __FILE__, __LINE__ );
    safe_fread( v->br.v, sizeof(real), bytesForVolume, model, __FILE__, __LINE__ );
    safe_fread( v->br.w, sizeof(real), bytesForVolume, model, __FILE__, __LINE__ );
    
    /* initialize density array  */
    safe_fread(      rho, sizeof(real), bytesForVolume, model, __FILE__, __LINE__ );

    /* close model file */
    safe_fclose ( "velocitymodel.bin", model, __FILE__, __LINE__ );     

    /* 
     * The elements on the HALO are not initialized here because they're going to be
     * exchanged at the beginning of the wave propagator 
     */

#else

    /* initalize velocity components */
    memset( v->tl.u, 0, cellsInArray );
    memset( v->tl.v, 0, cellsInArray );
    memset( v->tl.w, 0, cellsInArray );
    memset( v->tr.u, 0, cellsInArray );
    memset( v->tr.v, 0, cellsInArray );
    memset( v->tr.w, 0, cellsInArray );
    memset( v->bl.u, 0, cellsInArray );
    memset( v->bl.v, 0, cellsInArray );
    memset( v->bl.w, 0, cellsInArray );
    memset( v->br.u, 0, cellsInArray );
    memset( v->br.v, 0, cellsInArray );
    memset( v->br.w, 0, cellsInArray );
 
#endif // end of DDO_NOT_PERFOM_IO clause

    log_info("Initial velocity model loaded successfully");
};

void write_snapshot(char *folder,
                    int suffix,
                    real *data,
                    const integer numberOfCells,
                    const int nfields)
{
    char fname[300];
    sprintf(fname,"%s/snapshot.%05d.bin", folder, suffix);
    
    FILE *snapshot = safe_fopen(fname,"wb", __FILE__, __LINE__ );
    
    safe_fwrite(data, sizeof(real), numberOfCells * nfields, snapshot, __FILE__, __LINE__ );
    
    if ( fclose(snapshot)!=0)
        log_error ("can not close file %s", fname);
                /* not a critical error ? */
};

void read_snapshot(char *folder,
                   int suffix,
                   real *data,
                   const integer numberOfCells,
                   const int nfields)
{
    char fname[300];
    sprintf(fname,"%s/snapshot.%05d.bin", folder, suffix);
    
    FILE *snapshot = safe_fopen(fname,"rb", __FILE__, __LINE__ );
    
    safe_fread(data, sizeof(real), numberOfCells * nfields, snapshot, __FILE__, __LINE__ );
    
    if ( fclose(snapshot)!=0 || unlink(fname) != 0)
    {
        log_error ("Error closing or unliking file %s", fname);
                abort();
    }
};

void propagate_shot ( time_d        direction,
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
    /* local variables */
    integer cellsInHalo;   // number of cells being exchanged using MPI 
    int     rank;          // mpi local rank
    
    
    /* Initialize local variables */
    cellsInHalo = 10;  // (dimmz + 2*HALO) * (dimmx + 2*HALO) * HALO * sizeof(real); 
    MPI_Comm_rank ( MPI_COMM_WORLD, &rank );


    for(int t=0; t < timesteps; t++)
    {
        // if ( t%stacki == 0 && direction == BACKWARD) read_snapshot(folder, ntbwd-t, dataflush, datalen);
        
        // Phase 1 computation for velocities
        velocity_propagator(v, s, coeffs, rho, dt, dzi, dxi, dyi, nz0, nzf, nx0, nxf, ny0, ny0+HALO, dimmz, dimmx);
        velocity_propagator(v, s, coeffs, rho, dt, dzi, dxi, dyi, nz0, nzf, nx0, nxf, nyf-HALO, nyf, dimmz, dimmx);
        
        // Boundary exchange
        exchange_velocity_boundaries( &v, cellsInHalo, nyf, ny0);
        
        // Phase 2 computation for velocities
        velocity_propagator(v, s, coeffs, rho, dt, dzi, dxi, dyi, nz0, nzf, nx0, nxf, ny0+HALO, nyf-HALO, dimmz, dimmx);
        
        // Phase 1 computation for stresses
        stress_propagator(s, v, coeffs, rho, dt, dzi, dxi, dyi, nz0, nzf, nx0, nxf, ny0, ny0+HALO, dimmz, dimmx);
        stress_propagator(s, v, coeffs, rho, dt, dzi, dxi, dyi, nz0, nzf, nx0, nxf, nyf-HALO, nyf, dimmz, dimmx);
        
        // Boundary exchange
        exchange_stress_boundaries( &s, cellsInHalo, nyf, ny0);
        
        // P2 computation
        stress_propagator(s, v, coeffs, rho, dt, dzi, dxi, dyi, nz0, nzf, nx0, nxf, ny0+HALO, nyf-HALO, dimmz, dimmx);      
        
        // if ( t%stacki == 0 && direction == FORWARD) write_snapshot(folder, ntbwd-t, dataflush, datalen);
    }
    
    log_info ( "Kernel finished successfully.");
};

/*
NAME:exchange_boundaries
PURPOSE: data exchanges between the boundary layers of the analyzed volume

v                   (in) struct containing velocity arrays (4 points / cell x 3 components / point = 12 arrays)
numElement          (in) Number of elements to exchange
idxt                (in) identifier related to the folder
nyf                 (in) final plane to be exchanged
ny0                 (in) intial plane to be exchanged

RETURN none
*/
void exchange_velocity_boundaries ( v_t *v, 
                                    integer nCells, 
                                    integer nyf, 
                                    integer ny0 )
{
    log_info ("Exchanging velocity boundaries. Number of cells " I "", nCells);

    const uint64_t UNUSED(idxA) = (nyf-HALO)* nCells;
    const uint64_t UNUSED(idxB) = nyf       * nCells;
    const uint64_t UNUSED(idxC) = ny0       * nCells;

    int rank, ranksize;
    MPI_Comm_rank ( MPI_COMM_WORLD, &rank    );
    MPI_Comm_size ( MPI_COMM_WORLD, &ranksize);

    log_info ("         Rank %d out of rank %d", rank, ranksize );

    if ( rank != 0 )  //task to exchange velocities boundaries
    {
        EXCHANGE( &v->tl.u[0], &v->tl.u[0], rank-1, rank, nCells );
    
        /*
        //v tl u v w
        exchange_buffer ( &v->tl.u[0], &v->tl.u[0], rank-1, numberOfCells);
        exchange_buffer ( &v->tl.v[0], &v->tl.v[0], rank-1, numberOfCells);
        exchange_buffer ( &v->tl.w[0], &v->tl.w[0], rank-1, numberOfCells);
        
        //v tr u v w
        exchange_buffer ( &v->tr.u[0], &v->tr.u[0], rank-1, numberOfCells);
        exchange_buffer ( &v->tr.v[0], &v->tr.v[0], rank-1, numberOfCells);
        exchange_buffer ( &v->tr.w[0], &v->tr.w[0], rank-1, numberOfCells);
        
        //v b u v w
        exchange_buffer ( &v->bl.u[0], &v->bl.u[0], rank-1, numberOfCells);
        exchange_buffer ( &v->bl.v[0], &v->bl.v[0], rank-1, numberOfCells);
        exchange_buffer ( &v->bl.w[0], &v->bl.w[0], rank-1, numberOfCells);
        
        //v b1 u v w
        exchange_buffer ( &v->br.u[0], &v->br.u[0], rank-1, numberOfCells);
        exchange_buffer ( &v->br.v[0], &v->br.v[0], rank-1, numberOfCells);
        exchange_buffer ( &v->br.w[0], &v->br.w[0], rank-1, numberOfCells);
        */
    }

    if ( rank != ranksize -1 )  //task to exchange stress boundaries
    {
        EXCHANGE( &v->tl.u[0], &v->tl.u[0], rank+1, rank, nCells );
  
        /*
        //v tl u v w
        exchange_buffer ( &v->tl.u[0], &v->tl.u[0], rank+1, numberOfCells);
        exchange_buffer ( &v->tl.v[0], &v->tl.v[0], rank+1, numberOfCells);
        exchange_buffer ( &v->tl.w[0], &v->tl.w[0], rank+1, numberOfCells);
        
        //v tr u v w
        exchange_buffer ( &v->tr.u[0], &v->tr.u[0], rank+1, numberOfCells);
        exchange_buffer ( &v->tr.v[0], &v->tr.v[0], rank+1, numberOfCells);
        exchange_buffer ( &v->tr.w[0], &v->tr.w[0], rank+1, numberOfCells);
        
        //v b u v w
        exchange_buffer ( &v->bl.u[0], &v->bl.u[0], rank+1, numberOfCells);
        exchange_buffer ( &v->bl.v[0], &v->bl.v[0], rank+1, numberOfCells);
        exchange_buffer ( &v->bl.w[0], &v->bl.w[0], rank+1, numberOfCells);
        
        //v b1 u v w
        exchange_buffer ( &v->br.u[0], &v->br.u[0], rank+1, numberOfCells);
        exchange_buffer ( &v->br.v[0], &v->br.v[0], rank+1, numberOfCells);
        exchange_buffer ( &v->br.w[0], &v->br.w[0], rank+1, numberOfCells);
        */
    }

    log_info ("Velocity boundaries exchanged successfully");
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
void exchange_stress_boundaries ( s_t *s, 
                                  integer nCells, 
                                  integer nyf, 
                                  integer ny0 )
{
    log_info ("Exchanging stress boundaries");
 
    const uint64_t UNUSED(idxA) = (nyf-HALO)* nCells;
    const uint64_t UNUSED(idxB) = nyf       * nCells;
    const uint64_t UNUSED(idxC) = ny0       * nCells;

    int rank, ranksize;
    MPI_Comm_rank ( MPI_COMM_WORLD, &rank    );
    MPI_Comm_size ( MPI_COMM_WORLD, &ranksize);

    log_info ("         Rank %d out of rank %d", rank, ranksize );
  
    if ( rank != 0 ) //task to exchange velocities boundaries
    {
        EXCHANGE( &s->tl.zz[0], &s->tl.zz[0], rank-1, rank, nCells );

        /* 
        //s tl zz xz yz xx xy yy
        exchange_buffer ( &s->tl.zz[0], &s->tl.zz[0], rank-1, numberOfCells);
        exchange_buffer ( &s->tl.xz[0], &s->tl.xz[0], rank-1, numberOfCells);
        exchange_buffer ( &s->tl.yz[0], &s->tl.yz[0], rank-1, numberOfCells);
        exchange_buffer ( &s->tl.xx[0], &s->tl.xx[0], rank-1, numberOfCells);
        exchange_buffer ( &s->tl.xy[0], &s->tl.xy[0], rank-1, numberOfCells);
        exchange_buffer ( &s->tl.yy[0], &s->tl.yy[0], rank-1, numberOfCells);
        
        //s tr zz xz yz xx xy yy
        exchange_buffer ( &s->tr.zz[0], &s->tr.zz[0], rank-1, numberOfCells);
        exchange_buffer ( &s->tr.xz[0], &s->tr.xz[0], rank-1, numberOfCells);
        exchange_buffer ( &s->tr.yz[0], &s->tr.yz[0], rank-1, numberOfCells);
        exchange_buffer ( &s->tr.xx[0], &s->tr.xx[0], rank-1, numberOfCells);
        exchange_buffer ( &s->tr.xy[0], &s->tr.xy[0], rank-1, numberOfCells);
        exchange_buffer ( &s->tr.yy[0], &s->tr.yy[0], rank-1, numberOfCells);
        
        //s b zz xz yz xx xy yy
        exchange_buffer ( &s->bl.zz[0], &s->bl.zz[0], rank-1, numberOfCells);
        exchange_buffer ( &s->bl.xz[0], &s->bl.xz[0], rank-1, numberOfCells);
        exchange_buffer ( &s->bl.yz[0], &s->bl.yz[0], rank-1, numberOfCells);
        exchange_buffer ( &s->bl.xx[0], &s->bl.xx[0], rank-1, numberOfCells);
        exchange_buffer ( &s->bl.xy[0], &s->bl.xy[0], rank-1, numberOfCells);
        exchange_buffer ( &s->bl.yy[0], &s->bl.yy[0], rank-1, numberOfCells);
        
        //s b1 zz xz yz xx xy yy
        exchange_buffer ( &s->br.zz[0], &s->br.zz[0], rank-1, numberOfCells);
        exchange_buffer ( &s->br.xz[0], &s->br.xz[0], rank-1, numberOfCells);
        exchange_buffer ( &s->br.yz[0], &s->br.yz[0], rank-1, numberOfCells);
        exchange_buffer ( &s->br.xx[0], &s->br.xx[0], rank-1, numberOfCells);
        exchange_buffer ( &s->br.xy[0], &s->br.xy[0], rank-1, numberOfCells);
        exchange_buffer ( &s->br.yy[0], &s->br.yy[0], rank-1, numberOfCells);
        */
    }

    if ( rank != ranksize -1 )  //task to exchange stress boundaries
    {
        EXCHANGE( &s->tl.zz[0], &s->tl.zz[0], rank+1, rank, nCells );

        /* 
        //s tl zz xz yz xx xy yy
        exchange_buffer ( &s->tl.zz[0], &s->tl.zz[0], rank+1, numberOfCells);
        exchange_buffer ( &s->tl.xz[0], &s->tl.xz[0], rank+1, numberOfCells);
        exchange_buffer ( &s->tl.yz[0], &s->tl.yz[0], rank+1, numberOfCells);
        exchange_buffer ( &s->tl.xx[0], &s->tl.xx[0], rank+1, numberOfCells);
        exchange_buffer ( &s->tl.xy[0], &s->tl.xy[0], rank+1, numberOfCells);
        exchange_buffer ( &s->tl.yy[0], &s->tl.yy[0], rank+1, numberOfCells);
        
        //s tr zz xz yz xx xy yy
        exchange_buffer ( &s->tr.zz[0], &s->tr.zz[0], rank+1, numberOfCells);
        exchange_buffer ( &s->tr.xz[0], &s->tr.xz[0], rank+1, numberOfCells);
        exchange_buffer ( &s->tr.yz[0], &s->tr.yz[0], rank+1, numberOfCells);
        exchange_buffer ( &s->tr.xx[0], &s->tr.xx[0], rank+1, numberOfCells);
        exchange_buffer ( &s->tr.xy[0], &s->tr.xy[0], rank+1, numberOfCells);
        exchange_buffer ( &s->tr.yy[0], &s->tr.yy[0], rank+1, numberOfCells);
        
        //s b zz xz yz xx xy yy  
        exchange_buffer ( &s->bl.zz[0], &s->bl.zz[0], rank+1, numberOfCells);
        exchange_buffer ( &s->bl.xz[0], &s->bl.xz[0], rank+1, numberOfCells);
        exchange_buffer ( &s->bl.yz[0], &s->bl.yz[0], rank+1, numberOfCells);
        exchange_buffer ( &s->bl.xx[0], &s->bl.xx[0], rank+1, numberOfCells);
        exchange_buffer ( &s->bl.xy[0], &s->bl.xy[0], rank+1, numberOfCells);
        exchange_buffer ( &s->bl.yy[0], &s->bl.yy[0], rank+1, numberOfCells);
        
        //s b1 zz xz yz xx xy yy
        exchange_buffer ( &s->br.zz[0], &s->br.zz[0], rank+1, numberOfCells);
        exchange_buffer ( &s->br.xz[0], &s->br.xz[0], rank+1, numberOfCells);
        exchange_buffer ( &s->br.yz[0], &s->br.yz[0], rank+1, numberOfCells);
        exchange_buffer ( &s->br.xx[0], &s->br.xx[0], rank+1, numberOfCells);
        exchange_buffer ( &s->br.xy[0], &s->br.xy[0], rank+1, numberOfCells);
        exchange_buffer ( &s->br.yy[0], &s->br.yy[0], rank+1, numberOfCells);
        */
    }
    
    log_info ("Stress boundaries exchanged successfully");
};


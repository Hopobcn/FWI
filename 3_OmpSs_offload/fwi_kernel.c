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

 #ifdef DO_NOT_PERFOM_IO /* initalize velocity components */

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
		 char* fwipath = read_env_variable("FWIDIR");
		 int pathlen = strlen(fwipath) + 200;
     char modelname[pathlen];
     sprintf( modelname, "%s/InputModels/velocitymodel_%.2f.bin", fwipath, waveletFreq );

     fprintf(stderr, "Loading input model %s from disk (this could take a while)\n", modelname);

     FILE* model = safe_fopen( modelname, "rb", __FILE__, __LINE__ );

			double tstart=0.0, tend=0.0, iospeed=0.0;
			tstart = dtime();
    
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

		tend = dtime() - tstart;
		iospeed = ((numberOfCells * sizeof(real) * 12.f) / (1024.f * 1024.f)) / tend;

		fprintf(stderr, "Load initial model time: %lf seconds (%lf MiB/s)\n", tend, iospeed);
  
		/* close model file */
     safe_fclose ( "velocitymodel.bin", model, __FILE__, __LINE__ );

 #endif /* end of DDO_NOT_PERFOM_IO clause */
 };

/*
 * This funtion is used to store both the local preconditioner and gradient
 * fields during the kernel execution.
 */
void store_field (char *shotdir,
                  const int shotid,
									field_t type,
                  v_t *v,
                  const integer numberOfCells)
{
#ifdef DO_NOT_PERFOM_IO
  fprintf(stderr, "Warning: We are not doing any IO here (%s)\n", __FUNCTION__);

#else
	/* create field name */
	char fname[300];

	switch ( type )
	{
		case ( GRADIENT ):
		{
			fprintf(stderr, "Storing local gradient field\n");
			sprintf( fname, "%s/gradient_%05d.dat", shotdir, shotid );
			break;
		}
		case ( PRECONDITIONER ):
		{	
			fprintf(stderr, "Storing local preconditioner field\n");
			sprintf( fname, "%s/precond_%05d.dat" , shotdir, shotid );
			break;
		}
		default:
		{	
			fprintf(stderr, "Invalid field type identificator\n");
			abort();
		}
	}

  FILE *ffile = safe_fopen(fname,"wb", __FILE__, __LINE__ );

	double tstart=0.0, tend=0.0, iospeed=0.0;
	tstart = dtime();
	
	safe_fwrite( v->tr.u, sizeof(real), numberOfCells, ffile, __FILE__, __LINE__ );
	safe_fwrite( v->tr.v, sizeof(real), numberOfCells, ffile, __FILE__, __LINE__ );
	safe_fwrite( v->tr.w, sizeof(real), numberOfCells, ffile, __FILE__, __LINE__ );
	
	safe_fwrite( v->tl.u, sizeof(real), numberOfCells, ffile, __FILE__, __LINE__ );
	safe_fwrite( v->tl.v, sizeof(real), numberOfCells, ffile, __FILE__, __LINE__ );
	safe_fwrite( v->tl.w, sizeof(real), numberOfCells, ffile, __FILE__, __LINE__ );

	safe_fwrite( v->br.u, sizeof(real), numberOfCells, ffile, __FILE__, __LINE__ );
	safe_fwrite( v->br.v, sizeof(real), numberOfCells, ffile, __FILE__, __LINE__ );
	safe_fwrite( v->br.w, sizeof(real), numberOfCells, ffile, __FILE__, __LINE__ );
	
	safe_fwrite( v->bl.u, sizeof(real), numberOfCells, ffile, __FILE__, __LINE__ );
	safe_fwrite( v->bl.v, sizeof(real), numberOfCells, ffile, __FILE__, __LINE__ );
	safe_fwrite( v->bl.w, sizeof(real), numberOfCells, ffile, __FILE__, __LINE__ );

	tend = dtime() - tstart;
	iospeed = ((numberOfCells * sizeof(real) * 12.f) / (1024.f * 1024.f)) / tend;

	fprintf(stderr, "Write field time: %lf seconds (%lf MiB/s)\n", tend, iospeed);
  
	if ( fclose(ffile)!=0)
      fprintf(stderr,"Error closing file %s\n", fname);

#endif
};

/*
 * Saves the complete velocity field to disk.
 */
void write_snapshot(char *folder,
                    int suffix,
                    v_t *v,
                    const integer numberOfCells)
{
#ifdef DO_NOT_PERFOM_IO
  fprintf(stderr, "Warning: We are not doing any IO here (%s)\n", __FUNCTION__);

#else
	fprintf(stderr, "Writing snapshot\n");
  char fname[300];
  sprintf(fname,"%s/snapshot.%05d.bin", folder, suffix);

  FILE *snapshot = safe_fopen(fname,"wb", __FILE__, __LINE__ );

	double tstart=0.0, tend=0.0, iospeed=0.0;
	tstart = dtime();

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

	tend = dtime() - tstart;
	iospeed = ((numberOfCells * sizeof(real) * 12.f) / (1024.f * 1024.f)) / tend;

	fprintf(stderr, "Write snapshot time: %lf seconds (%lf MiB/s)\n", tend, iospeed);

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
#ifdef DO_NOT_PERFOM_IO
  fprintf(stderr, "Warning: We are not doing any IO here (%s)\n", __FUNCTION__);

#else
	fprintf(stderr, "Reading snapshot\n");
  char fname[300];
  sprintf(fname,"%s/snapshot.%05d.bin", folder, suffix);

  FILE *snapshot = safe_fopen(fname,"rb", __FILE__, __LINE__ );

	double tstart=0.0, tend=0.0, iospeed=0.0;
	tstart = dtime();
	
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

	tend = dtime() - tstart;
	iospeed = ((numberOfCells * sizeof(real) * 12.f) / (1024.f * 1024.f)) / tend;
  
	fprintf(stderr, "Read snapshot time: %lf seconds (%lf MiB/s)\n", tend, iospeed);
	
	if ( fclose(snapshot)!=0 )
      fprintf(stderr,"Error closing file %s\n", fname);

#endif
};


void propagate_shot ( time_d       direction,
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

	/* locate myself */
	int rank, ranksize;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &ranksize);
	fprintf(stderr, "propagate shot:  rank %d out of %d\n", rank, ranksize );



    for(int t=0; t < timesteps; t++)
    {
      fprintf(stderr, "Computing %d-th timestep\n", t);

      /* perform IO */
      if ( t%stacki == 0 && direction == BACKWARD) read_snapshot(folder, ntbwd-t, &v, datalen);

			/* ------------------------------------------------------------------------------ */
			/*                      VELOCITY COMPUTATION                                      */
			/* ------------------------------------------------------------------------------ */
      
			/* Phase 1. Computation of the left-most planes of the domain */
      velocity_propagator(v, s, coeffs, rho, dt, dzi, dxi, dyi,
                          nz0 +   HALO,
                          nzf -   HALO,
                          nx0 +   HALO,
                          nxf -   HALO,
                          ny0 +   HALO,
                          ny0 + 2*HALO,
                          dimmz, dimmx);

      /* Phase 1. Computation of the right-most planes of the domain */
      velocity_propagator(v, s, coeffs, rho, dt, dzi, dxi, dyi,
                          nz0 +   HALO,
                          nzf -   HALO,
                          nx0 +   HALO,
                          nxf -   HALO,
                          nyf - 2*HALO,
                          nyf -   HALO,
                          dimmz, dimmx);
      
			/* Boundary exchange for velocity values */
			if ( t < timesteps -1 )
      	int val = exchange_velocity_boundaries( &v, datalen, ny0, nyf, rank, ranksize);

      /* Phase 2. Computation of the central planes. */
      velocity_propagator(v, s, coeffs, rho, dt, dzi, dxi, dyi,
                          nz0 +   HALO,
                          nzf -   HALO,
                          nx0 +   HALO,
                          nxf -   HALO,
                          ny0 +   HALO,
                          nyf -   HALO,
                          dimmz, dimmx);


			/* ------------------------------------------------------------------------------ */
			/*                        STRESS COMPUTATION                                      */
			/* ------------------------------------------------------------------------------ */

      /* Phase 1. Computation of the left-most planes of the domain */
      stress_propagator ( s, v, coeffs, rho, dt, dzi, dxi, dyi, 
                          nz0 +   HALO,
                          nzf -   HALO,
                          nx0 +   HALO,
                          nxf -   HALO,
                          ny0 +   HALO,
                          ny0 + 2*HALO,
													dimmz, dimmx);
      
      /* Phase 1. Computation of the right-most planes of the domain */
			stress_propagator ( s, v, coeffs, rho, dt, dzi, dxi, dyi, 
                          nz0 +   HALO,
                          nzf -   HALO,
                          nx0 +   HALO,
                          nxf -   HALO,
                          nyf - 2*HALO,
                          nyf -   HALO,
													dimmz, dimmx);

      /* Boundary exchange for stress values */
      // exchange_stress_boundaries( &s, plane_size, rank, numTasks, nyf, ny0);
			if ( t < timesteps -1 )
      	int val = exchange_velocity_boundaries( &v, datalen, ny0, nyf, rank, ranksize);

      /* Phase 2 computation. Central planes of the domain */
			stress_propagator ( s, v, coeffs, rho, dt, dzi, dxi, dyi, 
                          nz0 +   HALO,
                          nzf -   HALO,
                          nx0 +   HALO,
                          nxf -   HALO,
                          ny0 +   HALO,
                          nyf -   HALO,
													dimmz, dimmx);
      	
			/* perform IO */
      if ( t%stacki == 0 && direction == FORWARD) write_snapshot(folder, ntbwd-t, &v, datalen);
    }
};

/*
 * dest is my local rank
 */
void exchange_buffer( real* A, real* B, int sendcount, int dest )
{
	MPI_Status status;
	int        tag    = 100;
	int res;
	
	void* sendbuf = (void *) A;
	void* recvbuf = (void *) B;

	res = MPI_Sendrecv( sendbuf, sendcount, MPI_FLOAT, dest, tag, 
								      recvbuf, sendcount, MPI_FLOAT, dest, tag,
								      MPI_COMM_WORLD, &status );

	if ( res != MPI_SUCCESS )
		fprintf(stderr, "%s: error\n", __FUNCTION__ );
};


int exchange_velocity_boundaries( v_t *v, 
																	const integer numberOfCells,
																	const integer ny0, 
																	const integer nyf,
																	int rank, 
																	int ranksize)
{

	int returnCode = 0;
  const integer message_size = 10;
	const integer rbsp = 1000; /* right boundary starting point */

  if ( ranksize == 1)
    return returnCode;

  if (rank != 0)
    // task to exchange velocities boundaries
    // #pragma omp task label(exchange vel boundaries backward) tied (con esta funciona)
    #pragma omp task in(v->tl.u[rbsp;message_size]) out(v->tl.u[0;message_size]) label(exchange vel boundaries)
    {
      //v tl u v w
			exchange_buffer( &v->tl.u[rbsp], &v->tl.u[0], message_size, rank-1);
			exchange_buffer( &v->tl.v[rbsp], &v->tl.v[0], message_size, rank-1);
      exchange_buffer( &v->tl.w[rbsp], &v->tl.w[0], message_size, rank-1);
      //v tr u v w
      exchange_buffer( &v->tr.u[rbsp], &v->tr.u[0], message_size, rank-1);
      exchange_buffer( &v->tr.v[rbsp], &v->tr.v[0], message_size, rank-1);
      exchange_buffer( &v->tr.w[rbsp], &v->tr.w[0], message_size, rank-1);
      //v b u v w
      exchange_buffer( &v->bl.u[rbsp], &v->bl.u[0], message_size, rank-1);
      exchange_buffer( &v->bl.v[rbsp], &v->bl.v[0], message_size, rank-1);
      exchange_buffer( &v->bl.w[rbsp], &v->bl.w[0], message_size, rank-1);
      //v b1 u v w
      exchange_buffer( &v->br.u[rbsp], &v->br.u[0], message_size, rank-1);
      exchange_buffer( &v->br.v[rbsp], &v->br.v[0], message_size, rank-1);
      exchange_buffer( &v->br.w[rbsp], &v->br.w[0], message_size, rank-1);
    }

  if (rank != ranksize -1)
    // task to exchange stress boundaries
    #pragma omp task in(v->tl.u[rbsp;message_size]) out(v->tl.u[0;message_size]) label(exchange stress boundaries)
    {
      //v tl u v w
      exchange_buffer( &v->tl.u[rbsp], &v->tl.u[0], message_size, rank+1);
      exchange_buffer( &v->tl.v[rbsp], &v->tl.v[0], message_size, rank+1);
      exchange_buffer( &v->tl.w[rbsp], &v->tl.w[0], message_size, rank+1);
      //v tr u v w
      exchange_buffer( &v->tr.u[rbsp], &v->tr.u[0], message_size, rank+1);
      exchange_buffer( &v->tr.v[rbsp], &v->tr.v[0], message_size, rank+1);
      exchange_buffer( &v->tr.w[rbsp], &v->tr.w[0], message_size, rank+1);
      //v b u v w
      exchange_buffer( &v->bl.u[rbsp], &v->bl.u[0], message_size, rank+1);
      exchange_buffer( &v->bl.v[rbsp], &v->bl.v[0], message_size, rank+1);
      exchange_buffer( &v->bl.w[rbsp], &v->bl.w[0], message_size, rank+1);
      //v b1 u v w
      exchange_buffer( &v->br.u[rbsp], &v->br.u[0], message_size, rank+1);
      exchange_buffer( &v->br.v[rbsp], &v->br.v[0], message_size, rank+1);
      exchange_buffer( &v->br.w[rbsp], &v->br.w[0], message_size, rank+1);
    }

  #pragma omp taskwait

  return returnCode;
}


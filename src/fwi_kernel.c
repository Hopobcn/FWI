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

/*
 * Initializes an array of length "length" to a random number.
 */
void set_array_to_random_real( real* restrict array, const integer length)
{
    const real randvalue = rand() / (1.0 * RAND_MAX);

    print_debug("Array is being initialized to %f", randvalue);

#if defined(_OPENACC)
    #pragma acc kernels copyin(array[0:length])
#endif
    for( integer i = 0; i < length; i++ )
        array[i] = randvalue;
}

/*
 * Initializes an array of length "length" to a constant floating point value.
 */
void set_array_to_constant( real* restrict array, const real value, const integer length)
{
#if defined(_OPENACC)
    #pragma acc kernels copyin(array[0:length])
#endif
    for( integer i = 0; i < length; i++ )
        array[i] = value;
}

void check_memory_shot( const dim_t dim,
                        coeff_t *c,
                        s_t     *s,
                        v_t     *v,
                        real    *rho)
{
#if defined(DEBUG)
    print_debug("Checking memory shot values");

    int numberOfCells = dim.pitch * dim.xsize * dim.ysize;

    real UNUSED(value);
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

        value = rho[i];
    }
#endif /* end of pragma DEBUG */
};


void alloc_memory_shot( const extent_t req,
                        dim_t   *dim,
                        coeff_t *c,
                        s_t     *s,
                        v_t     *v,
                        real    **rho)
{
    PUSH_RANGE

    //const integer size = numberOfCells * sizeof(real);

    //print_debug("ptr size = " I " bytes ("I" elements)", size, numberOfCells);

    /* allocate coefficients */
    c->c11 = (real*) malloc3d_host(dim, ALIGN_REAL, HALO, req);
    c->c12 = (real*) malloc3d_host(dim, ALIGN_REAL, HALO, req);
    c->c13 = (real*) malloc3d_host(dim, ALIGN_REAL, HALO, req);
    c->c14 = (real*) malloc3d_host(dim, ALIGN_REAL, HALO, req);
    c->c15 = (real*) malloc3d_host(dim, ALIGN_REAL, HALO, req);
    c->c16 = (real*) malloc3d_host(dim, ALIGN_REAL, HALO, req);

    c->c22 = (real*) malloc3d_host(dim, ALIGN_REAL, HALO, req);
    c->c23 = (real*) malloc3d_host(dim, ALIGN_REAL, HALO, req);
    c->c24 = (real*) malloc3d_host(dim, ALIGN_REAL, HALO, req);
    c->c25 = (real*) malloc3d_host(dim, ALIGN_REAL, HALO, req);
    c->c26 = (real*) malloc3d_host(dim, ALIGN_REAL, HALO, req);

    c->c33 = (real*) malloc3d_host(dim, ALIGN_REAL, HALO, req);
    c->c34 = (real*) malloc3d_host(dim, ALIGN_REAL, HALO, req);
    c->c35 = (real*) malloc3d_host(dim, ALIGN_REAL, HALO, req);
    c->c36 = (real*) malloc3d_host(dim, ALIGN_REAL, HALO, req);

    c->c44 = (real*) malloc3d_host(dim, ALIGN_REAL, HALO, req);
    c->c45 = (real*) malloc3d_host(dim, ALIGN_REAL, HALO, req);
    c->c46 = (real*) malloc3d_host(dim, ALIGN_REAL, HALO, req);

    c->c55 = (real*) malloc3d_host(dim, ALIGN_REAL, HALO, req);
    c->c56 = (real*) malloc3d_host(dim, ALIGN_REAL, HALO, req);
    c->c66 = (real*) malloc3d_host(dim, ALIGN_REAL, HALO, req);

    /* allocate velocity components */
    v->tl.u = (real*) malloc3d_host(dim, ALIGN_REAL, HALO, req);
    v->tl.v = (real*) malloc3d_host(dim, ALIGN_REAL, HALO, req);
    v->tl.w = (real*) malloc3d_host(dim, ALIGN_REAL, HALO, req);

    v->tr.u = (real*) malloc3d_host(dim, ALIGN_REAL, HALO, req);
    v->tr.v = (real*) malloc3d_host(dim, ALIGN_REAL, HALO, req);
    v->tr.w = (real*) malloc3d_host(dim, ALIGN_REAL, HALO, req);

    v->bl.u = (real*) malloc3d_host(dim, ALIGN_REAL, HALO, req);
    v->bl.v = (real*) malloc3d_host(dim, ALIGN_REAL, HALO, req);
    v->bl.w = (real*) malloc3d_host(dim, ALIGN_REAL, HALO, req);

    v->br.u = (real*) malloc3d_host(dim, ALIGN_REAL, HALO, req);
    v->br.v = (real*) malloc3d_host(dim, ALIGN_REAL, HALO, req);
    v->br.w = (real*) malloc3d_host(dim, ALIGN_REAL, HALO, req);

    /* allocate stress components   */
    s->tl.zz = (real*) malloc3d_host(dim, ALIGN_REAL, HALO, req);
    s->tl.xz = (real*) malloc3d_host(dim, ALIGN_REAL, HALO, req);
    s->tl.yz = (real*) malloc3d_host(dim, ALIGN_REAL, HALO, req);
    s->tl.xx = (real*) malloc3d_host(dim, ALIGN_REAL, HALO, req);
    s->tl.xy = (real*) malloc3d_host(dim, ALIGN_REAL, HALO, req);
    s->tl.yy = (real*) malloc3d_host(dim, ALIGN_REAL, HALO, req);

    s->tr.zz = (real*) malloc3d_host(dim, ALIGN_REAL, HALO, req);
    s->tr.xz = (real*) malloc3d_host(dim, ALIGN_REAL, HALO, req);
    s->tr.yz = (real*) malloc3d_host(dim, ALIGN_REAL, HALO, req);
    s->tr.xx = (real*) malloc3d_host(dim, ALIGN_REAL, HALO, req);
    s->tr.xy = (real*) malloc3d_host(dim, ALIGN_REAL, HALO, req);
    s->tr.yy = (real*) malloc3d_host(dim, ALIGN_REAL, HALO, req);

    s->bl.zz = (real*) malloc3d_host(dim, ALIGN_REAL, HALO, req);
    s->bl.xz = (real*) malloc3d_host(dim, ALIGN_REAL, HALO, req);
    s->bl.yz = (real*) malloc3d_host(dim, ALIGN_REAL, HALO, req);
    s->bl.xx = (real*) malloc3d_host(dim, ALIGN_REAL, HALO, req);
    s->bl.xy = (real*) malloc3d_host(dim, ALIGN_REAL, HALO, req);
    s->bl.yy = (real*) malloc3d_host(dim, ALIGN_REAL, HALO, req);

    s->br.zz = (real*) malloc3d_host(dim, ALIGN_REAL, HALO, req);
    s->br.xz = (real*) malloc3d_host(dim, ALIGN_REAL, HALO, req);
    s->br.yz = (real*) malloc3d_host(dim, ALIGN_REAL, HALO, req);
    s->br.xx = (real*) malloc3d_host(dim, ALIGN_REAL, HALO, req);
    s->br.xy = (real*) malloc3d_host(dim, ALIGN_REAL, HALO, req);
    s->br.yy = (real*) malloc3d_host(dim, ALIGN_REAL, HALO, req);

    /* allocate density array       */
    *rho = (real*) malloc3d_host(dim, ALIGN_REAL, HALO, req);

#if defined(_OPENACC)

    const integer datalen = dim->pitch * dim->xsize * dim->ysize;

    const real* rrho  = *rho;
    
    coeff_t c_h = *c;
    coeff_t c_d;

    c_d.c11 = (real*) malloc3d_device(dim, ALIGN_REAL, HALO, req, (void*)c->c11);
    c_d.c12 = (real*) malloc3d_device(dim, ALIGN_REAL, HALO, req, (void*)c->c12);
    c_d.c13 = (real*) malloc3d_device(dim, ALIGN_REAL, HALO, req, (void*)c->c13);
    c_d.c14 = (real*) malloc3d_device(dim, ALIGN_REAL, HALO, req, (void*)c->c14);
    c_d.c15 = (real*) malloc3d_device(dim, ALIGN_REAL, HALO, req, (void*)c->c15);
    c_d.c16 = (real*) malloc3d_device(dim, ALIGN_REAL, HALO, req, (void*)c->c16);

    c_d.c22 = (real*) malloc3d_device(dim, ALIGN_REAL, HALO, req, (void*)c->c22);
    c_d.c23 = (real*) malloc3d_device(dim, ALIGN_REAL, HALO, req, (void*)c->c23);
    c_d.c24 = (real*) malloc3d_device(dim, ALIGN_REAL, HALO, req, (void*)c->c24);
    c_d.c25 = (real*) malloc3d_device(dim, ALIGN_REAL, HALO, req, (void*)c->c25);
    c_d.c26 = (real*) malloc3d_device(dim, ALIGN_REAL, HALO, req, (void*)c->c26);

    c_d.c33 = (real*) malloc3d_device(dim, ALIGN_REAL, HALO, req, (void*)c->c33);
    c_d.c34 = (real*) malloc3d_device(dim, ALIGN_REAL, HALO, req, (void*)c->c34);
    c_d.c35 = (real*) malloc3d_device(dim, ALIGN_REAL, HALO, req, (void*)c->c35);
    c_d.c36 = (real*) malloc3d_device(dim, ALIGN_REAL, HALO, req, (void*)c->c36);

    c_d.c44 = (real*) malloc3d_device(dim, ALIGN_REAL, HALO, req, (void*)c->c44);
    c_d.c45 = (real*) malloc3d_device(dim, ALIGN_REAL, HALO, req, (void*)c->c45);
    c_d.c46 = (real*) malloc3d_device(dim, ALIGN_REAL, HALO, req, (void*)c->c46);

    c_d.c55 = (real*) malloc3d_device(dim, ALIGN_REAL, HALO, req, (void*)c->c55);
    c_d.c56 = (real*) malloc3d_device(dim, ALIGN_REAL, HALO, req, (void*)c->c56);
    c_d.c66 = (real*) malloc3d_device(dim, ALIGN_REAL, HALO, req, (void*)c->c66);
    
   //#pragma acc enter data create(cc)
   //#pragma acc enter data create(cc.c11[:datalen])
   //#pragma acc enter data create(cc.c12[:datalen])
   //#pragma acc enter data create(cc.c13[:datalen])
   //#pragma acc enter data create(cc.c14[:datalen])
   //#pragma acc enter data create(cc.c15[:datalen])
   //#pragma acc enter data create(cc.c16[:datalen])
   //#pragma acc enter data create(cc.c22[:datalen])
   //#pragma acc enter data create(cc.c23[:datalen])
   //#pragma acc enter data create(cc.c24[:datalen])
   //#pragma acc enter data create(cc.c25[:datalen])
   //#pragma acc enter data create(cc.c26[:datalen])
   //#pragma acc enter data create(cc.c33[:datalen])
   //#pragma acc enter data create(cc.c34[:datalen])
   //#pragma acc enter data create(cc.c35[:datalen])
   //#pragma acc enter data create(cc.c36[:datalen])
   //#pragma acc enter data create(cc.c44[:datalen])
   //#pragma acc enter data create(cc.c45[:datalen])
   //#pragma acc enter data create(cc.c46[:datalen])
   //#pragma acc enter data create(cc.c55[:datalen])
   //#pragma acc enter data create(cc.c56[:datalen])
   //#pragma acc enter data create(cc.c66[:datalen])

    v_t v_d;

    /* allocate velocity components */
    v_d.tl.u = (real*) malloc3d_device(dim, ALIGN_REAL, HALO, req, (void*)v->tl.u);
    v_d.tl.v = (real*) malloc3d_device(dim, ALIGN_REAL, HALO, req, (void*)v->tl.v);
    v_d.tl.w = (real*) malloc3d_device(dim, ALIGN_REAL, HALO, req, (void*)v->tl.w);
           
    v_d.tr.u = (real*) malloc3d_device(dim, ALIGN_REAL, HALO, req, (void*)v->tr.u);
    v_d.tr.v = (real*) malloc3d_device(dim, ALIGN_REAL, HALO, req, (void*)v->tr.v);
    v_d.tr.w = (real*) malloc3d_device(dim, ALIGN_REAL, HALO, req, (void*)v->tr.w);
           
    v_d.bl.u = (real*) malloc3d_device(dim, ALIGN_REAL, HALO, req, (void*)v->bl.u);
    v_d.bl.v = (real*) malloc3d_device(dim, ALIGN_REAL, HALO, req, (void*)v->bl.v);
    v_d.bl.w = (real*) malloc3d_device(dim, ALIGN_REAL, HALO, req, (void*)v->bl.w);
           
    v_d.br.u = (real*) malloc3d_device(dim, ALIGN_REAL, HALO, req, (void*)v->br.u);
    v_d.br.v = (real*) malloc3d_device(dim, ALIGN_REAL, HALO, req, (void*)v->br.v);
    v_d.br.w = (real*) malloc3d_device(dim, ALIGN_REAL, HALO, req, (void*)v->br.w);


   //v_t vv = *v;
   //#pragma acc enter data copyin(vv)
   //#pragma acc enter data create(vv.tl.u[:datalen])
   //#pragma acc enter data create(vv.tl.v[:datalen])
   //#pragma acc enter data create(vv.tl.w[:datalen])
   //#pragma acc enter data create(vv.tr.u[:datalen])
   //#pragma acc enter data create(vv.tr.v[:datalen])
   //#pragma acc enter data create(vv.tr.w[:datalen])
   //#pragma acc enter data create(vv.bl.u[:datalen])
   //#pragma acc enter data create(vv.bl.v[:datalen])
   //#pragma acc enter data create(vv.bl.w[:datalen])
   //#pragma acc enter data create(vv.br.u[:datalen])
   //#pragma acc enter data create(vv.br.v[:datalen])
   //#pragma acc enter data create(vv.br.w[:datalen])

    s_t s_d;
    /* allocate stress components   */
    s_d.tl.zz = (real*) malloc3d_device(dim, ALIGN_REAL, HALO, req, (void*) s->tl.zz);
    s_d.tl.xz = (real*) malloc3d_device(dim, ALIGN_REAL, HALO, req, (void*) s->tl.xz);
    s_d.tl.yz = (real*) malloc3d_device(dim, ALIGN_REAL, HALO, req, (void*) s->tl.yz);
    s_d.tl.xx = (real*) malloc3d_device(dim, ALIGN_REAL, HALO, req, (void*) s->tl.xx);
    s_d.tl.xy = (real*) malloc3d_device(dim, ALIGN_REAL, HALO, req, (void*) s->tl.xy);
    s_d.tl.yy = (real*) malloc3d_device(dim, ALIGN_REAL, HALO, req, (void*) s->tl.yy);
    
    s_d.tr.zz = (real*) malloc3d_device(dim, ALIGN_REAL, HALO, req, (void*) s->tr.zz);
    s_d.tr.xz = (real*) malloc3d_device(dim, ALIGN_REAL, HALO, req, (void*) s->tr.xz);
    s_d.tr.yz = (real*) malloc3d_device(dim, ALIGN_REAL, HALO, req, (void*) s->tr.yz);
    s_d.tr.xx = (real*) malloc3d_device(dim, ALIGN_REAL, HALO, req, (void*) s->tr.xx);
    s_d.tr.xy = (real*) malloc3d_device(dim, ALIGN_REAL, HALO, req, (void*) s->tr.xy);
    s_d.tr.yy = (real*) malloc3d_device(dim, ALIGN_REAL, HALO, req, (void*) s->tr.yy);
    
    s_d.bl.zz = (real*) malloc3d_device(dim, ALIGN_REAL, HALO, req, (void*) s->bl.zz);
    s_d.bl.xz = (real*) malloc3d_device(dim, ALIGN_REAL, HALO, req, (void*) s->bl.xz);
    s_d.bl.yz = (real*) malloc3d_device(dim, ALIGN_REAL, HALO, req, (void*) s->bl.yz);
    s_d.bl.xx = (real*) malloc3d_device(dim, ALIGN_REAL, HALO, req, (void*) s->bl.xx);
    s_d.bl.xy = (real*) malloc3d_device(dim, ALIGN_REAL, HALO, req, (void*) s->bl.xy);
    s_d.bl.yy = (real*) malloc3d_device(dim, ALIGN_REAL, HALO, req, (void*) s->bl.yy);
    
    s_d.br.zz = (real*) malloc3d_device(dim, ALIGN_REAL, HALO, req, (void*) s->br.zz);
    s_d.br.xz = (real*) malloc3d_device(dim, ALIGN_REAL, HALO, req, (void*) s->br.xz);
    s_d.br.yz = (real*) malloc3d_device(dim, ALIGN_REAL, HALO, req, (void*) s->br.yz);
    s_d.br.xx = (real*) malloc3d_device(dim, ALIGN_REAL, HALO, req, (void*) s->br.xx);
    s_d.br.xy = (real*) malloc3d_device(dim, ALIGN_REAL, HALO, req, (void*) s->br.xy);
    s_d.br.yy = (real*) malloc3d_device(dim, ALIGN_REAL, HALO, req, (void*) s->br.yy);


   //s_t ss = *s;
   //#pragma acc enter data copyin(ss)
   //#pragma acc enter data create(ss.tl.zz[:datalen])
   //#pragma acc enter data create(ss.tl.xz[:datalen])
   //#pragma acc enter data create(ss.tl.yz[:datalen])
   //#pragma acc enter data create(ss.tl.xx[:datalen])
   //#pragma acc enter data create(ss.tl.xy[:datalen])
   //#pragma acc enter data create(ss.tl.yy[:datalen])
   //#pragma acc enter data create(ss.tr.zz[:datalen])
   //#pragma acc enter data create(ss.tr.xz[:datalen])
   //#pragma acc enter data create(ss.tr.yz[:datalen])
   //#pragma acc enter data create(ss.tr.xx[:datalen])
   //#pragma acc enter data create(ss.tr.xy[:datalen])
   //#pragma acc enter data create(ss.tr.yy[:datalen])
   //#pragma acc enter data create(ss.bl.zz[:datalen])
   //#pragma acc enter data create(ss.bl.xz[:datalen])
   //#pragma acc enter data create(ss.bl.yz[:datalen])
   //#pragma acc enter data create(ss.bl.xx[:datalen])
   //#pragma acc enter data create(ss.bl.xy[:datalen])
   //#pragma acc enter data create(ss.bl.yy[:datalen])
   //#pragma acc enter data create(ss.br.zz[:datalen])
   //#pragma acc enter data create(ss.br.xz[:datalen])
   //#pragma acc enter data create(ss.br.yz[:datalen])
   //#pragma acc enter data create(ss.br.xx[:datalen])
   //#pragma acc enter data create(ss.br.xy[:datalen])
   //#pragma acc enter data create(ss.br.yy[:datalen])

    malloc3d_device(dim, ALIGN_REAL, HALO, req, (void*) *rho);

   //#pragma acc enter data create(rrho[:datalen])

#endif /* end of pragma _OPENACC */
    POP_RANGE
};

void free_memory_shot( coeff_t *c,
                       s_t     *s,
                       v_t     *v,
                       real    **rho)
{
    PUSH_RANGE

#if defined(_OPENACC)
    #pragma acc wait

   //#pragma acc exit data delete(c->c11)
   //#pragma acc exit data delete(c->c12)
   //#pragma acc exit data delete(c->c13)
   //#pragma acc exit data delete(c->c14)
   //#pragma acc exit data delete(c->c15)
   //#pragma acc exit data delete(c->c16)
   //#pragma acc exit data delete(c->c22)
   //#pragma acc exit data delete(c->c23)
   //#pragma acc exit data delete(c->c24)
   //#pragma acc exit data delete(c->c25)
   //#pragma acc exit data delete(c->c26)
   //#pragma acc exit data delete(c->c33)
   //#pragma acc exit data delete(c->c34)
   //#pragma acc exit data delete(c->c35)
   //#pragma acc exit data delete(c->c36)
   //#pragma acc exit data delete(c->c44)
   //#pragma acc exit data delete(c->c45)
   //#pragma acc exit data delete(c->c46)
   //#pragma acc exit data delete(c->c55)
   //#pragma acc exit data delete(c->c56)
   //#pragma acc exit data delete(c->c66)
   //#pragma acc exit data delete(c)

   //#pragma acc exit data delete(v->tl.u)
   //#pragma acc exit data delete(v->tl.v)
   //#pragma acc exit data delete(v->tl.w)
   //#pragma acc exit data delete(v->tr.u)
   //#pragma acc exit data delete(v->tr.v)
   //#pragma acc exit data delete(v->tr.w)
   //#pragma acc exit data delete(v->bl.u)
   //#pragma acc exit data delete(v->bl.v)
   //#pragma acc exit data delete(v->bl.w)
   //#pragma acc exit data delete(v->br.u)
   //#pragma acc exit data delete(v->br.v)
   //#pragma acc exit data delete(v->br.w)


   //#pragma acc exit data delete(s->tl.zz)
   //#pragma acc exit data delete(s->tl.xz)
   //#pragma acc exit data delete(s->tl.yz)
   //#pragma acc exit data delete(s->tl.xx)
   //#pragma acc exit data delete(s->tl.xy)
   //#pragma acc exit data delete(s->tl.yy)
   //#pragma acc exit data delete(s->tr.zz)
   //#pragma acc exit data delete(s->tr.xz)
   //#pragma acc exit data delete(s->tr.yz)
   //#pragma acc exit data delete(s->tr.xx)
   //#pragma acc exit data delete(s->tr.xy)
   //#pragma acc exit data delete(s->tr.yy)
   //#pragma acc exit data delete(s->bl.zz)
   //#pragma acc exit data delete(s->bl.xz)
   //#pragma acc exit data delete(s->bl.yz)
   //#pragma acc exit data delete(s->bl.xx)
   //#pragma acc exit data delete(s->bl.xy)
   //#pragma acc exit data delete(s->bl.yy)
   //#pragma acc exit data delete(s->br.zz)
   //#pragma acc exit data delete(s->br.xz)
   //#pragma acc exit data delete(s->br.yz)
   //#pragma acc exit data delete(s->br.xx)
   //#pragma acc exit data delete(s->br.xy)
   //#pragma acc exit data delete(s->br.yy)
   //#pragma acc exit data delete(s)
   //
   //const real* rrho  = *rho;
   //#pragma acc exit data delete(rrho)

    /* deallocate coefficients */
    free3d_device( (void*) c->c11 );
    free3d_device( (void*) c->c12 );
    free3d_device( (void*) c->c13 );
    free3d_device( (void*) c->c14 );
    free3d_device( (void*) c->c15 );
    free3d_device( (void*) c->c16 );

    free3d_device( (void*) c->c22 );
    free3d_device( (void*) c->c23 );
    free3d_device( (void*) c->c24 );
    free3d_device( (void*) c->c25 );
    free3d_device( (void*) c->c26 );
    free3d_device( (void*) c->c33 );

    free3d_device( (void*) c->c34 );
    free3d_device( (void*) c->c35 );
    free3d_device( (void*) c->c36 );

    free3d_device( (void*) c->c44 );
    free3d_device( (void*) c->c45 );
    free3d_device( (void*) c->c46 );

    free3d_device( (void*) c->c55 );
    free3d_device( (void*) c->c56 );

    free3d_device( (void*) c->c66 );

    /* deallocate velocity components */
    free3d_device( (void*) v->tl.u );
    free3d_device( (void*) v->tl.v );
    free3d_device( (void*) v->tl.w );

    free3d_device( (void*) v->tr.u );
    free3d_device( (void*) v->tr.v );
    free3d_device( (void*) v->tr.w );

    free3d_device( (void*) v->bl.u );
    free3d_device( (void*) v->bl.v );
    free3d_device( (void*) v->bl.w );

    free3d_device( (void*) v->br.u );
    free3d_device( (void*) v->br.v );
    free3d_device( (void*) v->br.w );

    /* deallocate stres components   */
    free3d_device( (void*) s->tl.zz );
    free3d_device( (void*) s->tl.xz );
    free3d_device( (void*) s->tl.yz );
    free3d_device( (void*) s->tl.xx );
    free3d_device( (void*) s->tl.xy );
    free3d_device( (void*) s->tl.yy );

    free3d_device( (void*) s->tr.zz );
    free3d_device( (void*) s->tr.xz );
    free3d_device( (void*) s->tr.yz );
    free3d_device( (void*) s->tr.xx );
    free3d_device( (void*) s->tr.xy );
    free3d_device( (void*) s->tr.yy );

    free3d_device( (void*) s->bl.zz );
    free3d_device( (void*) s->bl.xz );
    free3d_device( (void*) s->bl.yz );
    free3d_device( (void*) s->bl.xx );
    free3d_device( (void*) s->bl.xy );
    free3d_device( (void*) s->bl.yy );

    free3d_device( (void*) s->br.zz );
    free3d_device( (void*) s->br.xz );
    free3d_device( (void*) s->br.yz );
    free3d_device( (void*) s->br.xx );
    free3d_device( (void*) s->br.xy );
    free3d_device( (void*) s->br.yy );

    /* deallocate density array       */
    free3d_device( (void*) *rho );

#endif /* end pragma _OPENACC */

    /* deallocate coefficients */
    free3d_host( (void*) c->c11 );
    free3d_host( (void*) c->c12 );
    free3d_host( (void*) c->c13 );
    free3d_host( (void*) c->c14 );
    free3d_host( (void*) c->c15 );
    free3d_host( (void*) c->c16 );

    free3d_host( (void*) c->c22 );
    free3d_host( (void*) c->c23 );
    free3d_host( (void*) c->c24 );
    free3d_host( (void*) c->c25 );
    free3d_host( (void*) c->c26 );
    free3d_host( (void*) c->c33 );

    free3d_host( (void*) c->c34 );
    free3d_host( (void*) c->c35 );
    free3d_host( (void*) c->c36 );

    free3d_host( (void*) c->c44 );
    free3d_host( (void*) c->c45 );
    free3d_host( (void*) c->c46 );

    free3d_host( (void*) c->c55 );
    free3d_host( (void*) c->c56 );

    free3d_host( (void*) c->c66 );

    /* deallocate velocity components */
    free3d_host( (void*) v->tl.u );
    free3d_host( (void*) v->tl.v );
    free3d_host( (void*) v->tl.w );

    free3d_host( (void*) v->tr.u );
    free3d_host( (void*) v->tr.v );
    free3d_host( (void*) v->tr.w );

    free3d_host( (void*) v->bl.u );
    free3d_host( (void*) v->bl.v );
    free3d_host( (void*) v->bl.w );

    free3d_host( (void*) v->br.u );
    free3d_host( (void*) v->br.v );
    free3d_host( (void*) v->br.w );

    /* deallocate stres components   */
    free3d_host( (void*) s->tl.zz );
    free3d_host( (void*) s->tl.xz );
    free3d_host( (void*) s->tl.yz );
    free3d_host( (void*) s->tl.xx );
    free3d_host( (void*) s->tl.xy );
    free3d_host( (void*) s->tl.yy );

    free3d_host( (void*) s->tr.zz );
    free3d_host( (void*) s->tr.xz );
    free3d_host( (void*) s->tr.yz );
    free3d_host( (void*) s->tr.xx );
    free3d_host( (void*) s->tr.xy );
    free3d_host( (void*) s->tr.yy );

    free3d_host( (void*) s->bl.zz );
    free3d_host( (void*) s->bl.xz );
    free3d_host( (void*) s->bl.yz );
    free3d_host( (void*) s->bl.xx );
    free3d_host( (void*) s->bl.xy );
    free3d_host( (void*) s->bl.yy );

    free3d_host( (void*) s->br.zz );
    free3d_host( (void*) s->br.xz );
    free3d_host( (void*) s->br.yz );
    free3d_host( (void*) s->br.xx );
    free3d_host( (void*) s->br.xy );
    free3d_host( (void*) s->br.yy );

    /* deallocate density array       */
    free3d_host( (void*) *rho );

    POP_RANGE
};

/*
 * Loads initial values from coeffs, stress and velocity.
 */
void load_initial_model ( const real    waveletFreq,
                          const dim_t dim,
                          coeff_t *c,
                          s_t     *s,
                          v_t     *v,
                          real    *rho)
{
    PUSH_RANGE

    const int numberOfCells = dim.pitch * dim.xsize * dim.ysize;

    /* initialize stress */
    set_array_to_constant( s->tl.zz, 0, numberOfCells);
    set_array_to_constant( s->tl.xz, 0, numberOfCells);
    set_array_to_constant( s->tl.yz, 0, numberOfCells);
    set_array_to_constant( s->tl.xx, 0, numberOfCells);
    set_array_to_constant( s->tl.xy, 0, numberOfCells);
    set_array_to_constant( s->tl.yy, 0, numberOfCells);
    set_array_to_constant( s->tr.zz, 0, numberOfCells);
    set_array_to_constant( s->tr.xz, 0, numberOfCells);
    set_array_to_constant( s->tr.yz, 0, numberOfCells);
    set_array_to_constant( s->tr.xx, 0, numberOfCells);
    set_array_to_constant( s->tr.xy, 0, numberOfCells);
    set_array_to_constant( s->tr.yy, 0, numberOfCells);
    set_array_to_constant( s->bl.zz, 0, numberOfCells);
    set_array_to_constant( s->bl.xz, 0, numberOfCells);
    set_array_to_constant( s->bl.yz, 0, numberOfCells);
    set_array_to_constant( s->bl.xx, 0, numberOfCells);
    set_array_to_constant( s->bl.xy, 0, numberOfCells);
    set_array_to_constant( s->bl.yy, 0, numberOfCells);
    set_array_to_constant( s->br.zz, 0, numberOfCells);
    set_array_to_constant( s->br.xz, 0, numberOfCells);
    set_array_to_constant( s->br.yz, 0, numberOfCells);
    set_array_to_constant( s->br.xx, 0, numberOfCells);
    set_array_to_constant( s->br.xy, 0, numberOfCells);
    set_array_to_constant( s->br.yy, 0, numberOfCells);

#if defined(DO_NOT_PERFORM_IO)

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

    /* initalize velocity components */
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

    /* initialize coefficients */
    set_array_to_constant( c->c11, 1.0, numberOfCells);
    set_array_to_constant( c->c12, 1.0, numberOfCells);
    set_array_to_constant( c->c13, 1.0, numberOfCells);
    set_array_to_constant( c->c14, 1.0, numberOfCells);
    set_array_to_constant( c->c15, 1.0, numberOfCells);
    set_array_to_constant( c->c16, 1.0, numberOfCells);
    set_array_to_constant( c->c22, 1.0, numberOfCells);
    set_array_to_constant( c->c23, 1.0, numberOfCells);
    set_array_to_constant( c->c24, 1.0, numberOfCells);
    set_array_to_constant( c->c25, 1.0, numberOfCells);
    set_array_to_constant( c->c26, 1.0, numberOfCells);
    set_array_to_constant( c->c33, 1.0, numberOfCells);
    set_array_to_constant( c->c34, 1.0, numberOfCells);
    set_array_to_constant( c->c35, 1.0, numberOfCells);
    set_array_to_constant( c->c36, 1.0, numberOfCells);
    set_array_to_constant( c->c44, 1.0, numberOfCells);
    set_array_to_constant( c->c45, 1.0, numberOfCells);
    set_array_to_constant( c->c46, 1.0, numberOfCells);
    set_array_to_constant( c->c55, 1.0, numberOfCells);
    set_array_to_constant( c->c56, 1.0, numberOfCells);
    set_array_to_constant( c->c66, 1.0, numberOfCells);

    /* initialize rho */
    set_array_to_constant( rho, 1.0, numberOfCells );

    /* local variables */
    double tstart_outer, tstart_inner;
    double tend_outer, tend_inner;
    double iospeed_inner, iospeed_outer;
    char modelname[300];

     /* open initial model, binary file */
    sprintf( modelname, "../data/inputmodels/velocitymodel_%.2f.bin", waveletFreq );

    print_info("Loading input model %s from disk (this could take a while)", modelname);

    /* start clock, take into account file opening */
    tstart_outer = dtime();
    FILE* model = safe_fopen( modelname, "rb", __FILE__, __LINE__ );

    /* start clock, do not take into account file opening */
    tstart_inner = dtime();

    int id;
#if defined(USE_MPI)
    MPI_Comm_rank( MPI_COMM_WORLD, &id );
#else
    id = 0;
#endif

    const integer bytesForVolume = numberOfCells * sizeof(real);

    /* seek to the correct position corresponding to id (0 or rank) */
    if (fseek ( model, bytesForVolume * id, SEEK_SET) != 0)
        print_error("fseek() failed to set the correct position");

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

    /* stop inner timer */
    tend_inner = dtime() - tstart_inner;

    /* stop timer and compute statistics */
    safe_fclose ( "velocitymodel.bin", model, __FILE__, __LINE__ );
    tend_outer = dtime() - tstart_outer;

    //fprintf(stderr, "Number of cells %d\n", numberOfCells);
    //fprintf(stderr, "sizeof real %lu\n", sizeof(real));
    //fprintf(stderr, "bytes %lf\n", numberOfCells * sizeof(real) * 12.f);

    iospeed_inner = ((numberOfCells * sizeof(real) * 12.f) / (1000.f * 1000.f)) / tend_inner;
    iospeed_outer = ((numberOfCells * sizeof(real) * 12.f) / (1000.f * 1000.f)) / tend_outer;

    print_stats("Initial velocity model loaded (%lf GB)", TOGB(numberOfCells * sizeof(real) * 12));
    print_stats("\tInner time %lf seconds (%lf MiB/s)", tend_inner, iospeed_inner);
    print_stats("\tOuter time %lf seconds (%lf MiB/s)", tend_outer, iospeed_outer);
    print_stats("\tDifference %lf seconds", tend_outer - tend_inner);

#if defined(_OPENACC)
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

    #pragma acc update device(vtlu[0:numberOfCells], vtlv[0:numberOfCells], vtlw[0:numberOfCells]) \
                       device(vtru[0:numberOfCells], vtrv[0:numberOfCells], vtrw[0:numberOfCells]) \
                       device(vblu[0:numberOfCells], vblv[0:numberOfCells], vblw[0:numberOfCells]) \
                       device(vbru[0:numberOfCells], vbrv[0:numberOfCells], vbrw[0:numberOfCells]) \
                       async(H2D)
#endif /* end of pragma _OPENACC */
#endif /* end of pragma DDO_NOT_PERFORM_IO clause */

    POP_RANGE
};


/*
 * Saves the complete velocity field to disk.
 */
void write_snapshot(char *folder,
                    int suffix,
                    v_t *v,
                    const dim_t dim)
{
    PUSH_RANGE

#if defined(DO_NOT_PERFORM_IO)
    print_info("We are not writing the snapshot here cause IO is not enabled!");
#else

    int domain, ndomains;
#if defined(USE_MPI)
    MPI_Comm_rank( MPI_COMM_WORLD, &domain );
    MPI_Comm_size( MPI_COMM_WORLD, &ndomains );
#else
    domain = 0; ndomains = 1;
#endif

    const integer cellsInVolume  = (dim.pitch) * (dim.xsize) * ( (dim.ysize-2*HALO)/ndomains );
    const integer cellsInHALOs   = (dim.pitch) * (dim.xsize) * (2*HALO);
    const integer numberOfCells  = cellsInVolume + cellsInHALOs;
    const integer bytesForVolume = cellsInVolume * sizeof(real);

#if defined(_OPENACC)
    #pragma acc update self(v->tr.u[0:numberOfCells], v->tr.v[0:numberOfCells], v->tr.w[0:numberOfCells]) \
                       self(v->tl.u[0:numberOfCells], v->tl.v[0:numberOfCells], v->tl.w[0:numberOfCells]) \
                       self(v->br.u[0:numberOfCells], v->br.v[0:numberOfCells], v->br.w[0:numberOfCells]) \
                       self(v->bl.u[0:numberOfCells], v->bl.v[0:numberOfCells], v->bl.w[0:numberOfCells])
#endif /* end pragma _OPENACC*/

    /* local variables */
    char fname[300];

    /* open snapshot file and write results */
    sprintf(fname,"%s/snapshot.%05d.bin", folder, suffix);

#if defined(LOG_IO_STATS)
    double tstart_outer = dtime();
#endif
    FILE *snapshot = safe_fopen(fname,"wb", __FILE__, __LINE__ );
#if defined(LOG_IO_STATS)
    double tstart_inner = dtime();
#endif

    /* seek to the correct position corresponding to domain(id) */
    if (fseek ( snapshot, bytesForVolume * domain, SEEK_SET) != 0)
        print_error("fseek() failed to set the correct position");

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

#if defined(LOG_IO_STATS)
    /* stop inner timer */
    double tend_inner = dtime();
#endif
    /* close file and stop outer timer */
    safe_fclose(fname, snapshot, __FILE__, __LINE__ );
#if defined(LOG_IO_STATS)
    double tend_outer = dtime();

    double iospeed_inner = (( (double) numberOfCells * sizeof(real) * 12.f) / (1000.f * 1000.f)) / (tend_inner - tstart_inner);
    double iospeed_outer = (( (double) numberOfCells * sizeof(real) * 12.f) / (1000.f * 1000.f)) / (tend_outer - tstart_outer);

    print_stats("Write snapshot (%lf GB)", TOGB(numberOfCells * sizeof(real) * 12));
    print_stats("\tInner time %lf seconds (%lf MB/s)", (tend_inner - tstart_inner), iospeed_inner);
    print_stats("\tOuter time %lf seconds (%lf MB/s)", (tend_outer - tstart_outer), iospeed_outer);
    print_stats("\tDifference %lf seconds", tend_outer - tend_inner);
#endif /* end pragma LOG_IO_STATS */
#endif /* end pragma DO_NOT_PERFORM_IO */

    POP_RANGE
};

/*
 * Reads the complete velocity field from disk.
 */
void read_snapshot(char *folder,
                   int suffix,
                   v_t *v,
                   const dim_t dim)
{
    PUSH_RANGE

#if defined(DO_NOT_PERFORM_IO)
    print_info("We are not reading the snapshot here cause IO is not enabled!");
#else
    /* local variables */
    char fname[300];

    /* open file and read snapshot */
    sprintf(fname,"%s/snapshot.%05d.bin", folder, suffix);

#if defined(LOG_IO_STATS)
    double tstart_outer = dtime();
#endif
    FILE *snapshot = safe_fopen(fname,"rb", __FILE__, __LINE__ );
#if defined(LOG_IO_STATS)
    double tstart_inner = dtime();
#endif

    int domain, ndomains;
#if defined(USE_MPI)
    MPI_Comm_rank( MPI_COMM_WORLD, &domain );
    MPI_Comm_size( MPI_COMM_WORLD, &ndomains );
#else
    domain = 0; ndomains = 1;
#endif

    const integer cellsInVolume  = (dim.pitch) * (dim.xsize) * ( (dim.ysize-2*HALO)/ndomains );
    const integer cellsInHALOs   = (dim.pitch) * (dim.xsize) * (2*HALO);
    const integer numberOfCells  = cellsInVolume + cellsInHALOs;
    const integer bytesForVolume = cellsInVolume * sizeof(real);

    /* seek to the correct position corresponding to rank */
    if (fseek ( snapshot, bytesForVolume * domain, SEEK_SET) != 0)
        print_error("fseek() failed to set the correct position");

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

#if defined(LOG_IO_STATS)
    /* stop inner timer */
    double tend_inner = dtime() - tstart_inner;
#endif
    /* close file and stop outer timer */
    safe_fclose(fname, snapshot, __FILE__, __LINE__ );
#if defined(LOG_IO_STATS)
    double tend_outer = dtime() - tstart_outer;

    double iospeed_inner = ((numberOfCells * sizeof(real) * 12.f) / (1000.f * 1000.f)) / tend_inner;
    double iospeed_outer = ((numberOfCells * sizeof(real) * 12.f) / (1000.f * 1000.f)) / tend_outer;

    print_stats("Read snapshot (%lf GB)", TOGB(numberOfCells * sizeof(real) * 12));
    print_stats("\tInner time %lf seconds (%lf MiB/s)", tend_inner, iospeed_inner);
    print_stats("\tOuter time %lf seconds (%lf MiB/s)", tend_outer, iospeed_outer);
    print_stats("\tDifference %lf seconds", tend_outer - tend_inner);
#endif

#if defined(_OPENACC)
    #pragma acc update device(v->tr.u[0:numberOfCells], v->tr.v[0:numberOfCells], v->tr.w[0:numberOfCells]) \
                       device(v->tl.u[0:numberOfCells], v->tl.v[0:numberOfCells], v->tl.w[0:numberOfCells]) \
                       device(v->br.u[0:numberOfCells], v->br.v[0:numberOfCells], v->br.w[0:numberOfCells]) \
                       device(v->bl.u[0:numberOfCells], v->bl.v[0:numberOfCells], v->bl.w[0:numberOfCells]) \
                       async(H2D)
#endif /* end pragma _OPENACC */
#endif /* end pragma DO_NOT_PERFORM_IO */

    POP_RANGE
};

void propagate_shot(time_d        direction,
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
                    real          *UNUSED(dataflush),
                    dim_t         dim)
{
    PUSH_RANGE

    double tglobal_start, tglobal_total = 0.0;
    double tstress_start, tstress_total = 0.0;
    double tvel_start, tvel_total = 0.0;
    double megacells = 0.0;

    for(int t=0; t < timesteps; t++)
    {
        PUSH_RANGE

        if( t % 10 == 0 ) print_info("Computing %d-th timestep", t);

        /* perform IO */
        if ( t%stacki == 0 && direction == BACKWARD) read_snapshot(folder, ntbwd-t, &v, dim);

        tglobal_start = dtime();

        /* wait read_snapshot H2D copies */
#if defined(_OPENACC)
        #pragma acc wait(H2D) if ( (t%stacki == 0 && direction == BACKWARD) || t==0 )
#endif

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
                            dim,
                            ONE_L);

        /* Phase 1. Computation of the right-most planes of the domain */
        velocity_propagator(v, s, coeffs, rho, dt, dzi, dxi, dyi,
                            nz0 +   HALO,
                            nzf -   HALO,
                            nx0 +   HALO,
                            nxf -   HALO,
                            nyf - 2*HALO,
                            nyf -   HALO,
                            dim,
                            ONE_R);

        /* Phase 2. Computation of the central planes. */
        tvel_start = dtime();

        velocity_propagator(v, s, coeffs, rho, dt, dzi, dxi, dyi,
                            nz0 +   HALO,
                            nzf -   HALO,
                            nx0 +   HALO,
                            nxf -   HALO,
                            ny0 +   HALO,
                            nyf -   HALO,
                            dim,
                            TWO);
#if defined(USE_MPI)
        const integer plane_size = dim.pitch * dim.xsize;
        /* Boundary exchange for velocity values */
        exchange_velocity_boundaries( v, plane_size, nyf, ny0);
#endif
#if defined(_OPENACC)
        #pragma acc wait(ONE_L, ONE_R, TWO)
#endif
        tvel_total += (dtime() - tvel_start);

        /* ------------------------------------------------------------------------------ */
        /*                        STRESS COMPUTATION                                      */
        /* ------------------------------------------------------------------------------ */

        /* Phase 1. Computation of the left-most planes of the domain */
        stress_propagator(s, v, coeffs, rho, dt, dzi, dxi, dyi,
                          nz0 +   HALO,
                          nzf -   HALO,
                          nx0 +   HALO,
                          nxf -   HALO,
                          ny0 +   HALO,
                          ny0 + 2*HALO,
                          dim,
                          ONE_L);

        /* Phase 1. Computation of the right-most planes of the domain */
        stress_propagator(s, v, coeffs, rho, dt, dzi, dxi, dyi,
                          nz0 +   HALO,
                          nzf -   HALO,
                          nx0 +   HALO,
                          nxf -   HALO,
                          nyf - 2*HALO,
                          nyf -   HALO,
                          dim,
                          ONE_R);

        /* Phase 2 computation. Central planes of the domain */
        tstress_start = dtime();

        stress_propagator(s, v, coeffs, rho, dt, dzi, dxi, dyi,
                          nz0 +   HALO,
                          nzf -   HALO,
                          nx0 +   HALO,
                          nxf -   HALO,
                          ny0 +   HALO,
                          nyf -   HALO,
                          dim,
                          TWO);

#if defined(USE_MPI)
        /* Boundary exchange for stress values */
        exchange_stress_boundaries( s, plane_size, nyf, ny0);
#endif

#if defined(_OPENACC)
        #pragma acc wait(ONE_L, ONE_R, TWO, H2D, D2H)
#endif
        tstress_total += (dtime() - tstress_start);

        tglobal_total += (dtime() - tglobal_start);

        /* perform IO */
        if ( t%stacki == 0 && direction == FORWARD) write_snapshot(folder, ntbwd-t, &v, dim);

#if defined(USE_MPI)
        MPI_Barrier( MPI_COMM_WORLD );
#endif
        POP_RANGE
    }

    /* compute some statistics */
    megacells = ((nzf - nz0) * (nxf - nx0) * (nyf - ny0)) / 1e6;
    tglobal_total /= (double) timesteps;
    tstress_total /= (double) timesteps;
    tvel_total    /= (double) timesteps;

    print_stats("Maingrid GLOBAL   computation took %lf seconds - %lf Mcells/s", tglobal_total, (2*megacells) / tglobal_total);
    print_stats("Maingrid STRESS   computation took %lf seconds - %lf Mcells/s", tstress_total,  megacells / tstress_total);
    print_stats("Maingrid VELOCITY computation took %lf seconds - %lf Mcells/s", tvel_total, megacells / tvel_total);

    POP_RANGE
};

#if defined(USE_MPI)
/*
NAME:exchange_boundaries
PURPOSE: data exchanges between the boundary layers of the analyzed volume

v                   (in) struct containing velocity arrays (4 points / cell x 3 components / point = 12 arrays)
plane_size          (in) Number of elements per plane to exchange
nyf                 (in) final plane to be exchanged
ny0                 (in) intial plane to be exchanged

RETURN none
*/
void exchange_velocity_boundaries ( v_t v,
                                    const integer plane_size,
                                    const integer nyf,
                                    const integer ny0 )
{
    PUSH_RANGE

    int     rank;          // mpi local rank
    int     nranks;        // num mpi ranks

    /* Initialize local variables */
    MPI_Comm_rank ( MPI_COMM_WORLD, &rank );
    MPI_Comm_size ( MPI_COMM_WORLD, &nranks );

    const integer num_planes = HALO;
    const integer nelems     = num_planes * plane_size;

    const integer left_recv  = ny0;
    const integer left_send  = ny0+HALO;

    const integer right_recv = nyf-HALO;
    const integer right_send = nyf-2*HALO;

    if ( rank != 0 )
    {
        // [RANK-1] <---> [RANK] communication
        EXCHANGE( &v.tl.u[left_send], &v.tl.u[left_recv], rank-1, rank, nelems );
        EXCHANGE( &v.tl.v[left_send], &v.tl.v[left_recv], rank-1, rank, nelems );
        EXCHANGE( &v.tl.w[left_send], &v.tl.w[left_recv], rank-1, rank, nelems );

        EXCHANGE( &v.tr.u[left_send], &v.tr.u[left_recv], rank-1, rank, nelems );
        EXCHANGE( &v.tr.v[left_send], &v.tr.v[left_recv], rank-1, rank, nelems );
        EXCHANGE( &v.tr.w[left_send], &v.tr.w[left_recv], rank-1, rank, nelems );

        EXCHANGE( &v.bl.u[left_send], &v.bl.u[left_recv], rank-1, rank, nelems );
        EXCHANGE( &v.bl.v[left_send], &v.bl.v[left_recv], rank-1, rank, nelems );
        EXCHANGE( &v.bl.w[left_send], &v.bl.w[left_recv], rank-1, rank, nelems );

        EXCHANGE( &v.br.u[left_send], &v.br.u[left_recv], rank-1, rank, nelems );
        EXCHANGE( &v.br.v[left_send], &v.br.v[left_recv], rank-1, rank, nelems );
        EXCHANGE( &v.br.w[left_send], &v.br.w[left_recv], rank-1, rank, nelems );
    }

    if ( rank != nranks -1 )  //task to exchange stress boundaries
    {
        //                [RANK] <---> [RANK+1] communication
        EXCHANGE( &v.tl.u[right_send], &v.tl.u[right_recv], rank+1, rank, nelems );
        EXCHANGE( &v.tl.v[right_send], &v.tl.v[right_recv], rank+1, rank, nelems );
        EXCHANGE( &v.tl.w[right_send], &v.tl.w[right_recv], rank+1, rank, nelems );

        EXCHANGE( &v.tr.u[right_send], &v.tr.u[right_recv], rank+1, rank, nelems );
        EXCHANGE( &v.tr.v[right_send], &v.tr.v[right_recv], rank+1, rank, nelems );
        EXCHANGE( &v.tr.w[right_send], &v.tr.w[right_recv], rank+1, rank, nelems );

        EXCHANGE( &v.bl.u[right_send], &v.bl.u[right_recv], rank+1, rank, nelems );
        EXCHANGE( &v.bl.v[right_send], &v.bl.v[right_recv], rank+1, rank, nelems );
        EXCHANGE( &v.bl.w[right_send], &v.bl.w[right_recv], rank+1, rank, nelems );

        EXCHANGE( &v.br.u[right_send], &v.br.u[right_recv], rank+1, rank, nelems );
        EXCHANGE( &v.br.v[right_send], &v.br.v[right_recv], rank+1, rank, nelems );
        EXCHANGE( &v.br.w[right_send], &v.br.w[right_recv], rank+1, rank, nelems );
    }

    POP_RANGE
};

/*
NAME:exchange_stress_boundaries
PURPOSE: data exchanges between the boundary layers of the analyzed volume

s                   (in) struct containing stress arrays (4 points / cell x 6 components / point = 24 arrays)
plane_size          (in) Number of elements per plane to exchange
rank                (in) rank id (CPU id)
nranks              (in) number of CPUs
nyf                 (in) final plane to be exchanged
ny0                 (in) intial plane to be exchanged

RETURN none
*/
void exchange_stress_boundaries ( s_t s,
                                  const integer plane_size,
                                  const integer nyf,
                                  const integer ny0 )
{
    PUSH_RANGE

    int     rank;          // mpi local rank
    int     nranks;        // num mpi ranks

    /* Initialize local variables */
    MPI_Comm_rank ( MPI_COMM_WORLD, &rank );
    MPI_Comm_size ( MPI_COMM_WORLD, &nranks );

    const integer num_planes = HALO;
    const integer nelems     = num_planes * plane_size;

    const integer left_recv  = ny0;
    const integer left_send  = ny0+HALO;

    const integer right_recv = nyf-HALO;
    const integer right_send = nyf-2*HALO;

    if ( rank != 0 )
    {
        // [RANK-1] <---> [RANK] communication
        EXCHANGE( &s.tl.zz[left_send], &s.tl.zz[left_recv], rank-1, rank, nelems );
        EXCHANGE( &s.tl.xz[left_send], &s.tl.xz[left_recv], rank-1, rank, nelems );
        EXCHANGE( &s.tl.yz[left_send], &s.tl.yz[left_recv], rank-1, rank, nelems );
        EXCHANGE( &s.tl.xx[left_send], &s.tl.xx[left_recv], rank-1, rank, nelems );
        EXCHANGE( &s.tl.xy[left_send], &s.tl.xy[left_recv], rank-1, rank, nelems );
        EXCHANGE( &s.tl.yy[left_send], &s.tl.yy[left_recv], rank-1, rank, nelems );

        EXCHANGE( &s.tr.zz[left_send], &s.tr.zz[left_recv], rank-1, rank, nelems );
        EXCHANGE( &s.tr.xz[left_send], &s.tr.xz[left_recv], rank-1, rank, nelems );
        EXCHANGE( &s.tr.yz[left_send], &s.tr.yz[left_recv], rank-1, rank, nelems );
        EXCHANGE( &s.tr.xx[left_send], &s.tr.xx[left_recv], rank-1, rank, nelems );
        EXCHANGE( &s.tr.xy[left_send], &s.tr.xy[left_recv], rank-1, rank, nelems );
        EXCHANGE( &s.tr.yy[left_send], &s.tr.yy[left_recv], rank-1, rank, nelems );

        EXCHANGE( &s.bl.zz[left_send], &s.bl.zz[left_recv], rank-1, rank, nelems );
        EXCHANGE( &s.bl.xz[left_send], &s.bl.xz[left_recv], rank-1, rank, nelems );
        EXCHANGE( &s.bl.yz[left_send], &s.bl.yz[left_recv], rank-1, rank, nelems );
        EXCHANGE( &s.bl.xx[left_send], &s.bl.xx[left_recv], rank-1, rank, nelems );
        EXCHANGE( &s.bl.xy[left_send], &s.bl.xy[left_recv], rank-1, rank, nelems );
        EXCHANGE( &s.bl.yy[left_send], &s.bl.yy[left_recv], rank-1, rank, nelems );

        EXCHANGE( &s.br.zz[left_send], &s.br.zz[left_recv], rank-1, rank, nelems );
        EXCHANGE( &s.br.xz[left_send], &s.br.xz[left_recv], rank-1, rank, nelems );
        EXCHANGE( &s.br.yz[left_send], &s.br.yz[left_recv], rank-1, rank, nelems );
        EXCHANGE( &s.br.xx[left_send], &s.br.xx[left_recv], rank-1, rank, nelems );
        EXCHANGE( &s.br.xy[left_send], &s.br.xy[left_recv], rank-1, rank, nelems );
        EXCHANGE( &s.br.yy[left_send], &s.br.yy[left_recv], rank-1, rank, nelems );
    }

    if ( rank != nranks-1 )
    {
        //                [RANK] <---> [RANK+1] communication
        EXCHANGE( &s.tl.zz[right_send], &s.tl.zz[right_recv], rank+1, rank, nelems );
        EXCHANGE( &s.tl.xz[right_send], &s.tl.xz[right_recv], rank+1, rank, nelems );
        EXCHANGE( &s.tl.yz[right_send], &s.tl.yz[right_recv], rank+1, rank, nelems );
        EXCHANGE( &s.tl.xx[right_send], &s.tl.xx[right_recv], rank+1, rank, nelems );
        EXCHANGE( &s.tl.xy[right_send], &s.tl.xy[right_recv], rank+1, rank, nelems );
        EXCHANGE( &s.tl.yy[right_send], &s.tl.yy[right_recv], rank+1, rank, nelems );

        EXCHANGE( &s.tr.zz[right_send], &s.tr.zz[right_recv], rank+1, rank, nelems );
        EXCHANGE( &s.tr.xz[right_send], &s.tr.xz[right_recv], rank+1, rank, nelems );
        EXCHANGE( &s.tr.yz[right_send], &s.tr.yz[right_recv], rank+1, rank, nelems );
        EXCHANGE( &s.tr.xx[right_send], &s.tr.xx[right_recv], rank+1, rank, nelems );
        EXCHANGE( &s.tr.xy[right_send], &s.tr.xy[right_recv], rank+1, rank, nelems );
        EXCHANGE( &s.tr.yy[right_send], &s.tr.yy[right_recv], rank+1, rank, nelems );

        EXCHANGE( &s.bl.zz[right_send], &s.bl.zz[right_recv], rank+1, rank, nelems );
        EXCHANGE( &s.bl.xz[right_send], &s.bl.xz[right_recv], rank+1, rank, nelems );
        EXCHANGE( &s.bl.yz[right_send], &s.bl.yz[right_recv], rank+1, rank, nelems );
        EXCHANGE( &s.bl.xx[right_send], &s.bl.xx[right_recv], rank+1, rank, nelems );
        EXCHANGE( &s.bl.xy[right_send], &s.bl.xy[right_recv], rank+1, rank, nelems );
        EXCHANGE( &s.bl.yy[right_send], &s.bl.yy[right_recv], rank+1, rank, nelems );

        EXCHANGE( &s.br.zz[right_send], &s.br.zz[right_recv], rank+1, rank, nelems );
        EXCHANGE( &s.br.xz[right_send], &s.br.xz[right_recv], rank+1, rank, nelems );
        EXCHANGE( &s.br.yz[right_send], &s.br.yz[right_recv], rank+1, rank, nelems );
        EXCHANGE( &s.br.xx[right_send], &s.br.xx[right_recv], rank+1, rank, nelems );
        EXCHANGE( &s.br.xy[right_send], &s.br.xy[right_recv], rank+1, rank, nelems );
        EXCHANGE( &s.br.yy[right_send], &s.br.yy[right_recv], rank+1, rank, nelems );
    }

    POP_RANGE
};
#endif /* end of pragma USE_MPI */


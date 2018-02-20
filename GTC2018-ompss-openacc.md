We include an experimental version of the new OmpSs/OpenACC interoperability.

The general idea is to use OmpSs to track the dependences between tasks and manage GPU memory & Streams while using OpenACC to generate GPU kernels without having to use CUDA.

The following example shows an `openacc` OmpSs task:
```bash
const integer size = dimmz * dimmx * dimmy;

#pragma omp target device(openacc) copy_deps
#pragma omp task in( [size]rho, [size]sxptr, [size]syptr, [size]szptr ) inout( [size]vptr ) label(vcell_TL)
#pragma acc kernels deviceptr(rho, sxptr, syptr, szptr, vptr)
#pragma acc loop independent
for(integer y=ny0; y < nyf; y++)
{
    #pragma acc loop independent
    for(integer x=nx0; x < nxf; x++)
    {
        #pragma acc loop independent
        for(integer z=nz0; z < nzf; z++)
        {
            const real lrho = RHO_BL( rho,z,x,y,dimmz,dimmx );

            const real stx = STENCILX( sxptr,SX,dxi,z,x,y,dimmz,dimmx );
            const real sty = STENCILY( syptr,SY,dyi,z,x,y,dimmz,dimmx );
            const real stz = STENCILZ( szptr,SZ,dzi,z,x,y,dimmz,dimmx );

            vptr[IDX(z,x,y,dimmz,dimmx)] += (stx  + sty  + stz) * dt * lrho;
        }
    }
}
```
OmpSs `target device(openacc)` tells the compiler that the following `omp task` is a task that have to be executed in a GPU.

#### Current Limitations:
> At the moment, the user has to manually specify all the GPU pointers to the OpenACC parallel region using the `deviceptr` clause.
That way the OpenACC will assume that you (in this case the OmpSs runtime) manages the GPU memory and won't attempt to allocate/copy/deallocate your GPU arrays.

> Also, the clause `async` is already provided by the `mercurium` compiler (in OmpSs, the streams are managed by the runtime). 
> So providing additional `async` clauses could brake compilation.

> To avoid coherence issues, we are limited to only One GPU (the user have to provide `NX_ARGS="--gpus=1"` in every execution)

#### Build the example:
```bash
cd FWI-sol-ompss-acc
make
NX_ARGS="--gpus=1" ./fwi data/fwi_params.txt data/fwi_frequencies.profile.txt
```

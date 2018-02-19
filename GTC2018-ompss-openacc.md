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
We only support the `kernels` and `loop` directive with clauses `deviceptr`. 
Any other combination is not supported.
Also, we are restricted to only one GPU.

To build this example:
```bash
cd FWI-sol-ompss-acc
make -i
NX_ARGS="--gpus=1" ./fwi data/fwi_params.txt data/fwi_frequencies.profile.txt
```

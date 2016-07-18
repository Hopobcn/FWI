#include <cuda_runtime.h>
#include <cstdio>

//extern "C" {
//    #include "fwi_propagator.h"
//}


#define C0 1.2f
#define C1 1.4f
#define C2 1.6f
#define C3 1.8f


__device__ inline
int IDX (const int z, 
         const int x, 
         const int y, 
         const int dimmz, 
         const int dimmx)
{
    return (y*dimmx*dimmz) + (x*dimmz) + (z);
};

__device__ inline
float stencil_Z (const int off,
                 const float* __restrict__ ptr,
                 const float    dzi,
                 const int z,
                 const int x,
                 const int y,
                 const int dimmz,
                 const int dimmx)
{
    return  ((C0 * ( ptr[IDX(z  +off,x,y,dimmz,dimmx)] - ptr[IDX(z-1+off,x,y,dimmz,dimmx)]) +
              C1 * ( ptr[IDX(z+1+off,x,y,dimmz,dimmx)] - ptr[IDX(z-2+off,x,y,dimmz,dimmx)]) +
              C2 * ( ptr[IDX(z+2+off,x,y,dimmz,dimmx)] - ptr[IDX(z-3+off,x,y,dimmz,dimmx)]) +
              C3 * ( ptr[IDX(z+3+off,x,y,dimmz,dimmx)] - ptr[IDX(z-4+off,x,y,dimmz,dimmx)])) * dzi );
};

__device__ inline
float stencil_X(const int off,
                const float* __restrict__ ptr,
                const float dxi,
                const int z,
                const int x,
                const int y,
                const int dimmz,
                const int dimmx)
{
    return ((C0 * ( ptr[IDX(z,x  +off,y,dimmz,dimmx)] - ptr[IDX(z,x-1+off,y,dimmz,dimmx)]) +
             C1 * ( ptr[IDX(z,x+1+off,y,dimmz,dimmx)] - ptr[IDX(z,x-2+off,y,dimmz,dimmx)]) +
             C2 * ( ptr[IDX(z,x+2+off,y,dimmz,dimmx)] - ptr[IDX(z,x-3+off,y,dimmz,dimmx)]) +
             C3 * ( ptr[IDX(z,x+3+off,y,dimmz,dimmx)] - ptr[IDX(z,x-4+off,y,dimmz,dimmx)])) * dxi );
};

__device__ inline
float stencil_Y(const int off,
                const float* __restrict__ ptr,
                const float dyi,
                const int z,
                const int x,
                const int y,
                const int dimmz,
                const int dimmx)
{
    return ((C0 * ( ptr[IDX(z,x,y  +off,dimmz,dimmx)] - ptr[IDX(z,x,y-1+off,dimmz,dimmx)]) +
             C1 * ( ptr[IDX(z,x,y+1+off,dimmz,dimmx)] - ptr[IDX(z,x,y-2+off,dimmz,dimmx)]) +
             C2 * ( ptr[IDX(z,x,y+2+off,dimmz,dimmx)] - ptr[IDX(z,x,y-3+off,dimmz,dimmx)]) +
             C3 * ( ptr[IDX(z,x,y+3+off,dimmz,dimmx)] - ptr[IDX(z,x,y-4+off,dimmz,dimmx)])) * dyi );
};

/* -------------------------------------------------------------------- */
/*                     KERNELS FOR VELOCITY                             */
/* -------------------------------------------------------------------- */

__device__ inline
float rho_BL (const float* __restrict__ rho,
              const int z,
              const int x,
              const int y,
              const int dimmz,
              const int dimmx)
{
    return (2.0f / (rho[IDX(z,x,y,dimmz,dimmx)] + rho[IDX(z+1,x,y,dimmz,dimmx)]));
};

__device__ inline
float rho_TR (const float* __restrict__ rho,
              const int z,
              const int x,
              const int y,
              const int dimmz,
              const int dimmx)
{
    return (2.0f / (rho[IDX(z,x,y,dimmz,dimmx)] + rho[IDX(z,x+1,y,dimmz,dimmx)]));
};

__device__ inline
float rho_BR (const float* __restrict__ rho,
              const int z,
              const int x,
              const int y,
              const int dimmz,
              const int dimmx)
{
    return ( 8.0f/ ( rho[IDX(z  ,x  ,y  ,dimmz,dimmx)] +
                     rho[IDX(z+1,x  ,y  ,dimmz,dimmx)] +
                     rho[IDX(z  ,x+1,y  ,dimmz,dimmx)] +
                     rho[IDX(z  ,x  ,y+1,dimmz,dimmx)] +
                     rho[IDX(z  ,x+1,y+1,dimmz,dimmx)] +
                     rho[IDX(z+1,x+1,y  ,dimmz,dimmx)] +
                     rho[IDX(z+1,x  ,y+1,dimmz,dimmx)] +
                     rho[IDX(z+1,x+1,y+1,dimmz,dimmx)]) );
};

__device__ inline
float rho_TL (const float* __restrict__ rho,
              const int z,
              const int x,
              const int y,
              const int dimmz,
              const int dimmx)
{
    return (2.0f / (rho[IDX(z,x,y,dimmz,dimmx)] + rho[IDX(z,x,y+1,dimmz,dimmx)]));
};

#ifdef OPTIMIZED
template <const int HALO = 4, 
          const int BLOCK_DIM_X = 16,
          const int BLOCK_DIM_Y = 16>
__global__
__launch_bounds__(128, 16) 
void compute_component_vcell_TL_cuda_k ( float* __restrict__ vptr,
                                   const float* __restrict__ szptr,
                                   const float* __restrict__ sxptr,
                                   const float* __restrict__ syptr,
                                   const float* __restrict__ rho,
                                   const float           dt,
                                   const float           dzi,
                                   const float           dxi,
                                   const float           dyi,
                                   const int             nz0,
                                   const int             nzf,
                                   const int             nx0,
                                   const int             nxf,
                                   const int             ny0,
                                   const int             nyf,
                                   const int             SZ,
                                   const int             SX,
                                   const int             SY,
                                   const int             dimmz,
                                   const int             dimmx)
{
    for(int z = blockIdx.x * blockDim.x + threadIdx.x + nz0; 
            z < nzf; 
            z += gridDim.x * blockDim.x)
    {
        for(int x = blockIdx.y * blockDim.y + threadIdx.y + nx0; 
                x < nxf; 
                x += gridDim.y * blockDim.y)
        {
            __shared__ float sx_smem[BLOCK_DIM_Y][BLOCK_DIM_X+2*HALO];
            float sy_front1, sy_front2, sy_front3, sy_front4;
            float sy_back1, sy_back2, sy_back3, sy_back4;
            float sy_current;

            sy_back3   = syptr[IDX(z,x,ny0-HALO+0+SY,dimmz,dimmx)];
            sy_back2   = syptr[IDX(z,x,ny0-HALO+1+SY,dimmz,dimmx)];
            sy_back1   = syptr[IDX(z,x,ny0-HALO+2+SY,dimmz,dimmx)];
            sy_current = syptr[IDX(z,x,ny0-HALO+3+SY,dimmz,dimmx)];
            sy_front1  = syptr[IDX(z,x,ny0-HALO+4+SY,dimmz,dimmx)];
            sy_front2  = syptr[IDX(z,x,ny0-HALO+5+SY,dimmz,dimmx)];
            sy_front3  = syptr[IDX(z,x,ny0-HALO+6+SY,dimmz,dimmx)];
            sy_front4  = syptr[IDX(z,x,ny0-HALO+7+SY,dimmz,dimmx)];

            float rho_current, rho_front1;
            rho_front1 = rho[IDX(z,x,ny0,dimmz,dimmx)];

            for(int y = ny0; 
                    y < nyf; 
                    y++)
            {
                /////////// register tiling-advance plane ////////////////
                sy_back4   = sy_back3;
                sy_back3   = sy_back2;
                sy_back2   = sy_back1;
                sy_back1   = sy_current;
                sy_current = sy_front1;
                sy_front1  = sy_front2;
                sy_front2  = sy_front3;
                sy_front3  = sy_front4;
                sy_front4  = syptr[IDX(z,x,y+SY+HALO,dimmz,dimmx)];
                ///////////////////////
                rho_current = rho_front1;
                rho_front1  = rho[IDX(z,x,y+1,dimmz,dimmx)];
                //////////////////////////////////////////////////////////

                const int tx = threadIdx.x+HALO;
                const int ty = threadIdx.y;
                ///////////// intra-block communication///////////////////
                sx_smem[ty][tx] = sxptr[IDX(z,x+SX,y,dimmz,dimmx)];
                //if (threadIdx.y < HALO)
                //{
                //    sx_smem[threadIdx.y                 ][tx] = sxptr[IDX(z,x+SX-HALO,       y,dimmz,dimmx)];
                //    sx_smem[threadIdx.y+BLOCK_DIM_Y+HALO][tx] = sxptr[IDX(z,x+SX+BLOCK_DIM_X,y,dimmz,dimmx)];
                //}
                if (threadIdx.x < HALO)
                {
                    sx_smem[ty][threadIdx.x                 ] = sxptr[IDX(z,x+SX-HALO,       y,dimmz,dimmx)];
                    sx_smem[ty][threadIdx.x+BLOCK_DIM_X+HALO] = sxptr[IDX(z,x+SX+BLOCK_DIM_X,y,dimmz,dimmx)];
                }
                /////////////////////////////////////////////////////////
                
                ///////////// intra-warp communication /////////////////
                float sz_current = szptr[IDX(z+SZ,x,y,dimmz,dimmx)];
                float sz_front3 = __shfl_up(sz_current, 3);
                float sz_front2 = __shfl_up(sz_current, 2);
                float sz_front1 = __shfl_up(sz_current, 1);
                float sz_back1  = __shfl_down(sz_current, 1);
                float sz_back2  = __shfl_down(sz_current, 2);
                float sz_back3  = __shfl_down(sz_current, 3);
                float sz_back4  = __shfl_down(sz_current, 4);
                ////////////////////////////////////////////////////////
                __syncthreads();

                const float lrho = (2.0f / (rho_current + rho_front1));

                const float stz = ((C0 * ( sz_current - sz_back1 ) +
                                    C1 * ( sz_front1  - sz_back2 ) +
                                    C2 * ( sz_front2  - sz_back3 ) +
                                    C3 * ( sz_front3  - sz_back4 )) * dzi );
                
                const float stx = ((C0 * ( sx_smem[ty][tx  ] - sx_smem[ty][tx-1] ) +
                                    C1 * ( sx_smem[ty][tx+1] - sx_smem[ty][tx-2] ) +
                                    C2 * ( sx_smem[ty][tx+2] - sx_smem[ty][tx-3] ) +
                                    C3 * ( sx_smem[ty][tx+3] - sx_smem[ty][tx-4] )) * dxi );

                const float sty = ((C0 * ( sy_current - sy_back1 ) +
                                    C1 * ( sy_front1  - sy_back2 ) +
                                    C2 * ( sy_front2  - sy_back3 ) +
                                    C3 * ( sy_front3  - sy_back4 )) * dyi );

                vptr[IDX(z,x,y,dimmz,dimmx)] += (stx  + sty  + stz) * dt * lrho;
            }
        }
    }
}
#else
__global__
void compute_component_vcell_TL_cuda_k ( float* __restrict__ vptr,
                                   const float* __restrict__ szptr,
                                   const float* __restrict__ sxptr,
                                   const float* __restrict__ syptr,
                                   const float* __restrict__ rho,
                                   const float           dt,
                                   const float           dzi,
                                   const float           dxi,
                                   const float           dyi,
                                   const int             nz0,
                                   const int             nzf,
                                   const int             nx0,
                                   const int             nxf,
                                   const int             ny0,
                                   const int             nyf,
                                   const int             SZ,
                                   const int             SX,
                                   const int             SY,
                                   const int             dimmz,
                                   const int             dimmx)
{
    for(int z = blockIdx.x * blockDim.x + threadIdx.x + nz0; 
            z < nzf; 
            z += gridDim.x * blockDim.x)
    {
        for(int x = blockIdx.y * blockDim.y + threadIdx.y + nx0; 
                x < nxf; 
                x += gridDim.y * blockDim.y)
        {
            for(int y = ny0; 
                    y < nyf; 
                    y++)
            {
                const float lrho = rho_TL(rho, z, x, y, dimmz, dimmx);
                
                const float stx  = stencil_X( SX, sxptr, dxi, z, x, y, dimmz, dimmx);
                const float sty  = stencil_Y( SY, syptr, dyi, z, x, y, dimmz, dimmx);
                const float stz  = stencil_Z( SZ, szptr, dzi, z, x, y, dimmz, dimmx);
                
                vptr[IDX(z,x,y,dimmz,dimmx)] += (stx  + sty  + stz) * dt * lrho;
            }
        }
    }
}
#endif

extern "C"
void compute_component_vcell_TL_cuda ( float* vptr,
                                 const float* szptr,
                                 const float* sxptr,
                                 const float* syptr,
                                 const float* rho,
                                 const float  dt,
                                 const float  dzi,
                                 const float  dxi,
                                 const float  dyi,
                                 const int    nz0,
                                 const int    nzf,
                                 const int    nx0,
                                 const int    nxf,
                                 const int    ny0,
                                 const int    nyf,
                                 const int    SZ,
                                 const int    SX,
                                 const int    SY,
                                 const int    dimmz,
                                 const int    dimmx,
                                 void*        stream)
{
    const int block_dim_x = 32;
    const int block_dim_y = 4;


    dim3 grid_dim( ((nzf-nz0) + block_dim_x-1)/block_dim_x,
                   ((nxf-nx0) + block_dim_y-1)/block_dim_y,
                      1 );
    dim3 block_dim(block_dim_x, block_dim_y, 1);

    cudaStream_t s = (cudaStream_t) stream;

#ifdef OPTIMIZED
    compute_component_vcell_TL_cuda_k<4,block_dim_x,block_dim_y><<<grid_dim, block_dim, 0, s>>>
        (vptr, szptr, sxptr, syptr, rho, dt, dzi, dxi, dyi, 
         nz0, nzf, nx0, nxf, ny0, nyf, SZ, SX, SY, dimmz, dimmx);
#else
    compute_component_vcell_TL_cuda_k<<<grid_dim, block_dim, 0, s>>>
        (vptr, szptr, sxptr, syptr, rho, dt, dzi, dxi, dyi, 
         nz0, nzf, nx0, nxf, ny0, nyf, SZ, SX, SY, dimmz, dimmx);
#endif
};


#ifdef OPTIMIZED
template <const int HALO = 4, 
          const int BLOCK_DIM_X = 16,
          const int BLOCK_DIM_Y = 16>
__global__
__launch_bounds__(128, 16) 
void compute_component_vcell_TR_cuda_k ( float* __restrict__ vptr,
                                   const float* __restrict__ szptr,
                                   const float* __restrict__ sxptr,
                                   const float* __restrict__ syptr,
                                   const float* __restrict__ rho,
                                   const float           dt,
                                   const float           dzi,
                                   const float           dxi,
                                   const float           dyi,
                                   const int             nz0,
                                   const int             nzf,
                                   const int             nx0,
                                   const int             nxf,
                                   const int             ny0,
                                   const int             nyf,
                                   const int             SZ,
                                   const int             SX,
                                   const int             SY,
                                   const int             dimmz,
                                   const int             dimmx)
{
    for(int z = blockIdx.x * blockDim.x + threadIdx.x + nz0; 
            z < nzf; 
            z += gridDim.x * blockDim.x)
    {
        for(int x = blockIdx.y * blockDim.y + threadIdx.y + nx0; 
                x < nxf; 
                x += gridDim.y * blockDim.y)
        {
            __shared__ float sx_smem[BLOCK_DIM_Y][BLOCK_DIM_X+2*HALO];
            float sy_front1, sy_front2, sy_front3, sy_front4;
            float sy_back1, sy_back2, sy_back3, sy_back4;
            float sy_current;

            sy_back3   = syptr[IDX(z,x,ny0-HALO+0+SY,dimmz,dimmx)];
            sy_back2   = syptr[IDX(z,x,ny0-HALO+1+SY,dimmz,dimmx)];
            sy_back1   = syptr[IDX(z,x,ny0-HALO+2+SY,dimmz,dimmx)];
            sy_current = syptr[IDX(z,x,ny0-HALO+3+SY,dimmz,dimmx)];
            sy_front1  = syptr[IDX(z,x,ny0-HALO+4+SY,dimmz,dimmx)];
            sy_front2  = syptr[IDX(z,x,ny0-HALO+5+SY,dimmz,dimmx)];
            sy_front3  = syptr[IDX(z,x,ny0-HALO+6+SY,dimmz,dimmx)];
            sy_front4  = syptr[IDX(z,x,ny0-HALO+7+SY,dimmz,dimmx)];

            __shared__ float rho_smem[BLOCK_DIM_Y+1][BLOCK_DIM_X];

            for(int y = ny0; 
                    y < nyf; 
                    y++)
            {
                /////////// register tiling-advance plane ////////////////
                sy_back4   = sy_back3;
                sy_back3   = sy_back2;
                sy_back2   = sy_back1;
                sy_back1   = sy_current;
                sy_current = sy_front1;
                sy_front1  = sy_front2;
                sy_front2  = sy_front3;
                sy_front3  = sy_front4;
                sy_front4  = syptr[IDX(z,x,y+SY+HALO,dimmz,dimmx)];
                //////////////////////////////////////////////////////////

                const int tx = threadIdx.x+HALO;
                const int ty = threadIdx.y;
                ///////////// intra-block communication///////////////////
                sx_smem[ty][tx] = sxptr[IDX(z,x+SX,y,dimmz,dimmx)];
                //if (threadIdx.y < HALO)
                //{
                //    sx_smem[threadIdx.y                 ][tx] = sxptr[IDX(z,x+SX-HALO,       y,dimmz,dimmx)];
                //    sx_smem[threadIdx.y+BLOCK_DIM_Y+HALO][tx] = sxptr[IDX(z,x+SX+BLOCK_DIM_X,y,dimmz,dimmx)];
                //}
                if (threadIdx.x < HALO)
                {
                    sx_smem[ty][threadIdx.x                 ] = sxptr[IDX(z,x+SX-HALO,       y,dimmz,dimmx)];
                    sx_smem[ty][threadIdx.x+BLOCK_DIM_X+HALO] = sxptr[IDX(z,x+SX+BLOCK_DIM_X,y,dimmz,dimmx)];
                }
                /////////////////////////////////////////////////////////
                rho_smem[ty][threadIdx.x] = rho[IDX(z,x,y,dimmz,dimmx)];
                if (threadIdx.y < 1)
                {
                    rho_smem[ty+BLOCK_DIM_Y+1][threadIdx.x] = rho[IDX(z,x+1,y,dimmz,dimmx)];
                }
                /////////////////////////////////////////////////////////
                ///////////// intra-warp communication /////////////////
                float sz_current = szptr[IDX(z+SZ,x,y,dimmz,dimmx)];
                float sz_front3 = __shfl_up(sz_current, 3);
                float sz_front2 = __shfl_up(sz_current, 2);
                float sz_front1 = __shfl_up(sz_current, 1);
                float sz_back1  = __shfl_down(sz_current, 1);
                float sz_back2  = __shfl_down(sz_current, 2);
                float sz_back3  = __shfl_down(sz_current, 3);
                float sz_back4  = __shfl_down(sz_current, 4);
                ////////////////////////////////////////////////////////
                __syncthreads();

                const float lrho = (2.0f / (rho_smem[ty][threadIdx.x] + rho_smem[ty+1][threadIdx.x]));

                const float stz = ((C0 * ( sz_current - sz_back1 ) +
                                    C1 * ( sz_front1  - sz_back2 ) +
                                    C2 * ( sz_front2  - sz_back3 ) +
                                    C3 * ( sz_front3  - sz_back4 )) * dzi );
                
                const float stx = ((C0 * ( sx_smem[ty][tx  ] - sx_smem[ty][tx-1] ) +
                                    C1 * ( sx_smem[ty][tx+1] - sx_smem[ty][tx-2] ) +
                                    C2 * ( sx_smem[ty][tx+2] - sx_smem[ty][tx-3] ) +
                                    C3 * ( sx_smem[ty][tx+3] - sx_smem[ty][tx-4] )) * dxi );

                const float sty = ((C0 * ( sy_current - sy_back1 ) +
                                    C1 * ( sy_front1  - sy_back2 ) +
                                    C2 * ( sy_front2  - sy_back3 ) +
                                    C3 * ( sy_front3  - sy_back4 )) * dyi );

                vptr[IDX(z,x,y,dimmz,dimmx)] += (stx  + sty  + stz) * dt * lrho;
            }
        }
    }
}
#else
__global__
void compute_component_vcell_TR_cuda_k ( float* __restrict__ vptr,
                                   const float* __restrict__ szptr,
                                   const float* __restrict__ sxptr,
                                   const float* __restrict__ syptr,
                                   const float* __restrict__ rho,
                                   const float           dt,
                                   const float           dzi,
                                   const float           dxi,
                                   const float           dyi,
                                   const int             nz0,
                                   const int             nzf,
                                   const int             nx0,
                                   const int             nxf,
                                   const int             ny0,
                                   const int             nyf,
                                   const int             SZ,
                                   const int             SX,
                                   const int             SY,
                                   const int             dimmz,
                                   const int             dimmx)
{
    for(int z = blockIdx.x * blockDim.x + threadIdx.x + nz0; 
            z < nzf; 
            z += gridDim.x * blockDim.x)
    {
        for(int x = blockIdx.y * blockDim.y + threadIdx.y + nx0; 
                x < nxf; 
                x += gridDim.y * blockDim.y)
        {
            for(int y = ny0; 
                    y < nyf; 
                    y++)
            {
                const float lrho = rho_TR(rho, z, x, y, dimmz, dimmx);
                
                const float stx  = stencil_X( SX, sxptr, dxi, z, x, y, dimmz, dimmx);
                const float sty  = stencil_Y( SY, syptr, dyi, z, x, y, dimmz, dimmx);
                const float stz  = stencil_Z( SZ, szptr, dzi, z, x, y, dimmz, dimmx);
                
                vptr[IDX(z,x,y,dimmz,dimmx)] += (stx  + sty  + stz) * dt * lrho;
            }
        }
    }
}
#endif

extern "C"
void compute_component_vcell_TR_cuda ( float* vptr,
                                 const float* szptr,
                                 const float* sxptr,
                                 const float* syptr,
                                 const float* rho,
                                 const float  dt,
                                 const float  dzi,
                                 const float  dxi,
                                 const float  dyi,
                                 const int    nz0,
                                 const int    nzf,
                                 const int    nx0,
                                 const int    nxf,
                                 const int    ny0,
                                 const int    nyf,
                                 const int    SZ,
                                 const int    SX,
                                 const int    SY,
                                 const int    dimmz,
                                 const int    dimmx,
                                 void*        stream)
{
    const int block_dim_x = 32;
    const int block_dim_y = 4;


    dim3 grid_dim( ((nzf-nz0) + block_dim_x-1)/block_dim_x,
                   ((nxf-nx0) + block_dim_y-1)/block_dim_y,
                      1 );
    dim3 block_dim(block_dim_x, block_dim_y, 1);

    cudaStream_t s = (cudaStream_t) stream;

#ifdef OPTIMIZED
    compute_component_vcell_TR_cuda_k<4,block_dim_x,block_dim_y><<<grid_dim, block_dim, 0, s>>>
        (vptr, szptr, sxptr, syptr, rho, dt, dzi, dxi, dyi, 
         nz0, nzf, nx0, nxf, ny0, nyf, SZ, SX, SY, dimmz, dimmx);
#else
    compute_component_vcell_TR_cuda_k<<<grid_dim, block_dim, 0, s>>>
        (vptr, szptr, sxptr, syptr, rho, dt, dzi, dxi, dyi, 
         nz0, nzf, nx0, nxf, ny0, nyf, SZ, SX, SY, dimmz, dimmx);
#endif
};


#ifdef OPTIMIZED
template <const int HALO = 4, 
          const int BLOCK_DIM_X = 16,
          const int BLOCK_DIM_Y = 16>
__global__
__launch_bounds__(128, 16) 
void compute_component_vcell_BR_cuda_k ( float* __restrict__ vptr,
                                   const float* __restrict__ szptr,
                                   const float* __restrict__ sxptr,
                                   const float* __restrict__ syptr,
                                   const float* __restrict__ rho,
                                   const float           dt,
                                   const float           dzi,
                                   const float           dxi,
                                   const float           dyi,
                                   const int             nz0,
                                   const int             nzf,
                                   const int             nx0,
                                   const int             nxf,
                                   const int             ny0,
                                   const int             nyf,
                                   const int             SZ,
                                   const int             SX,
                                   const int             SY,
                                   const int             dimmz,
                                   const int             dimmx)
{
    for(int z = blockIdx.x * blockDim.x + threadIdx.x + nz0; 
            z < nzf; 
            z += gridDim.x * blockDim.x)
    {
        for(int x = blockIdx.y * blockDim.y + threadIdx.y + nx0; 
                x < nxf; 
                x += gridDim.y * blockDim.y)
        {
            __shared__ float sx_smem[BLOCK_DIM_Y][BLOCK_DIM_X+2*HALO];
            float sy_front1, sy_front2, sy_front3, sy_front4;
            float sy_back1, sy_back2, sy_back3, sy_back4;
            float sy_current;

            sy_back3   = syptr[IDX(z,x,ny0-HALO+0+SY,dimmz,dimmx)];
            sy_back2   = syptr[IDX(z,x,ny0-HALO+1+SY,dimmz,dimmx)];
            sy_back1   = syptr[IDX(z,x,ny0-HALO+2+SY,dimmz,dimmx)];
            sy_current = syptr[IDX(z,x,ny0-HALO+3+SY,dimmz,dimmx)];
            sy_front1  = syptr[IDX(z,x,ny0-HALO+4+SY,dimmz,dimmx)];
            sy_front2  = syptr[IDX(z,x,ny0-HALO+5+SY,dimmz,dimmx)];
            sy_front3  = syptr[IDX(z,x,ny0-HALO+6+SY,dimmz,dimmx)];
            sy_front4  = syptr[IDX(z,x,ny0-HALO+7+SY,dimmz,dimmx)];

            __shared__ float rho_smem[BLOCK_DIM_Y+1][BLOCK_DIM_X+1];
            float rho_current, rho_front1;
            rho_front1 = rho[IDX(z,x,ny0,dimmz,dimmx)];

            for(int y = ny0; 
                    y < nyf; 
                    y++)
            {
                /////////// register tiling-advance plane ////////////////
                sy_back4   = sy_back3;
                sy_back3   = sy_back2;
                sy_back2   = sy_back1;
                sy_back1   = sy_current;
                sy_current = sy_front1;
                sy_front1  = sy_front2;
                sy_front2  = sy_front3;
                sy_front3  = sy_front4;
                sy_front4  = syptr[IDX(z,x,y+SY+HALO,dimmz,dimmx)];
                ///////////////////////
                rho_current = rho_front1;
                rho_front1  = rho[IDX(z,x,y+1,dimmz,dimmx)];
                //////////////////////////////////////////////////////////

                const int tx = threadIdx.x+HALO;
                const int ty = threadIdx.y;
                ///////////// intra-block communication///////////////////
                sx_smem[ty][tx] = sxptr[IDX(z,x+SX,y,dimmz,dimmx)];
                //if (threadIdx.y < HALO)
                //{
                //    sx_smem[threadIdx.y                 ][tx] = sxptr[IDX(z,x+SX-HALO,       y,dimmz,dimmx)];
                //    sx_smem[threadIdx.y+BLOCK_DIM_Y+HALO][tx] = sxptr[IDX(z,x+SX+BLOCK_DIM_X,y,dimmz,dimmx)];
                //}
                if (threadIdx.x < HALO)
                {
                    sx_smem[ty][threadIdx.x                 ] = sxptr[IDX(z,x+SX-HALO,       y,dimmz,dimmx)];
                    sx_smem[ty][threadIdx.x+BLOCK_DIM_X+HALO] = sxptr[IDX(z,x+SX+BLOCK_DIM_X,y,dimmz,dimmx)];
                }
                /////////////////////////////////////////////////////////
                rho_smem[ty][threadIdx.x] = rho[IDX(z,x,y,dimmz,dimmx)];
                if (threadIdx.y < 1)
                    rho_smem[ty+BLOCK_DIM_Y+1][threadIdx.x              ] = rho[IDX(z,x+BLOCK_DIM_Y+1,y,dimmz,dimmx)];
                if (threadIdx.x < 1)
                    rho_smem[ty              ][threadIdx.x+BLOCK_DIM_X+1] = rho[IDX(z+BLOCK_DIM_X+1,x,y,dimmz,dimmx)];
                /////////////////////////////////////////////////////////
                
                ///////////// intra-warp communication /////////////////
                float sz_current = szptr[IDX(z+SZ,x,y,dimmz,dimmx)];
                float sz_front3 = __shfl_up(sz_current, 3);
                float sz_front2 = __shfl_up(sz_current, 2);
                float sz_front1 = __shfl_up(sz_current, 1);
                float sz_back1  = __shfl_down(sz_current, 1);
                float sz_back2  = __shfl_down(sz_current, 2);
                float sz_back3  = __shfl_down(sz_current, 3);
                float sz_back4  = __shfl_down(sz_current, 4);
                ////////////////////////////////////////////////////////
                __syncthreads();

                const float lrho = (8.0f/ ( rho_current                       +
                                            rho_smem[ty  ][threadIdx.x+1]     +
                                            rho_smem[ty+1][threadIdx.x  ]     +
                                            rho_front1                        +
                                            rho[IDX(z  ,x+1,y+1,dimmz,dimmx)] +
                                            rho[IDX(z+1,x+1,y  ,dimmz,dimmx)] +
                                            rho[IDX(z+1,x  ,y+1,dimmz,dimmx)] +
                                            rho[IDX(z+1,x+1,y+1,dimmz,dimmx)]) );

                const float stz = ((C0 * ( sz_current - sz_back1 ) +
                                    C1 * ( sz_front1  - sz_back2 ) +
                                    C2 * ( sz_front2  - sz_back3 ) +
                                    C3 * ( sz_front3  - sz_back4 )) * dzi );
                
                const float stx = ((C0 * ( sx_smem[ty][tx  ] - sx_smem[ty][tx-1] ) +
                                    C1 * ( sx_smem[ty][tx+1] - sx_smem[ty][tx-2] ) +
                                    C2 * ( sx_smem[ty][tx+2] - sx_smem[ty][tx-3] ) +
                                    C3 * ( sx_smem[ty][tx+3] - sx_smem[ty][tx-4] )) * dxi );

                const float sty = ((C0 * ( sy_current - sy_back1 ) +
                                    C1 * ( sy_front1  - sy_back2 ) +
                                    C2 * ( sy_front2  - sy_back3 ) +
                                    C3 * ( sy_front3  - sy_back4 )) * dyi );

                vptr[IDX(z,x,y,dimmz,dimmx)] += (stx  + sty  + stz) * dt * lrho;
            }
        }
    }
}
#else
__global__
void compute_component_vcell_BR_cuda_k ( float* __restrict__ vptr,
                                   const float* __restrict__ szptr,
                                   const float* __restrict__ sxptr,
                                   const float* __restrict__ syptr,
                                   const float* __restrict__ rho,
                                   const float           dt,
                                   const float           dzi,
                                   const float           dxi,
                                   const float           dyi,
                                   const int             nz0,
                                   const int             nzf,
                                   const int             nx0,
                                   const int             nxf,
                                   const int             ny0,
                                   const int             nyf,
                                   const int             SZ,
                                   const int             SX,
                                   const int             SY,
                                   const int             dimmz,
                                   const int             dimmx)
{
    for(int z = blockIdx.x * blockDim.x + threadIdx.x + nz0; 
            z < nzf; 
            z += gridDim.x * blockDim.x)
    {
        for(int x = blockIdx.y * blockDim.y + threadIdx.y + nx0; 
                x < nxf; 
                x += gridDim.y * blockDim.y)
        {
            for(int y = ny0; 
                    y < nyf; 
                    y++)
            {
                const float lrho = rho_BR(rho, z, x, y, dimmz, dimmx);
                
                const float stx  = stencil_X( SX, sxptr, dxi, z, x, y, dimmz, dimmx);
                const float sty  = stencil_Y( SY, syptr, dyi, z, x, y, dimmz, dimmx);
                const float stz  = stencil_Z( SZ, szptr, dzi, z, x, y, dimmz, dimmx);
                
                vptr[IDX(z,x,y,dimmz,dimmx)] += (stx  + sty  + stz) * dt * lrho;
            }
        }
    }
}
#endif

extern "C"
void compute_component_vcell_BR_cuda ( float* vptr,
                                 const float* szptr,
                                 const float* sxptr,
                                 const float* syptr,
                                 const float* rho,
                                 const float  dt,
                                 const float  dzi,
                                 const float  dxi,
                                 const float  dyi,
                                 const int    nz0,
                                 const int    nzf,
                                 const int    nx0,
                                 const int    nxf,
                                 const int    ny0,
                                 const int    nyf,
                                 const int    SZ,
                                 const int    SX,
                                 const int    SY,
                                 const int    dimmz,
                                 const int    dimmx,
                                 void*        stream)
{
    const int block_dim_x = 32;
    const int block_dim_y = 4;


    dim3 grid_dim( ((nzf-nz0) + block_dim_x-1)/block_dim_x,
                   ((nxf-nx0) + block_dim_y-1)/block_dim_y,
                      1 );
    dim3 block_dim(block_dim_x, block_dim_y, 1);

    cudaStream_t s = (cudaStream_t) stream;

#ifdef OPTIMIZED
    compute_component_vcell_BR_cuda_k<4,block_dim_x,block_dim_y><<<grid_dim, block_dim, 0, s>>>
        (vptr, szptr, sxptr, syptr, rho, dt, dzi, dxi, dyi, 
         nz0, nzf, nx0, nxf, ny0, nyf, SZ, SX, SY, dimmz, dimmx);
#else
    compute_component_vcell_BR_cuda_k<<<grid_dim, block_dim, 0, s>>>
        (vptr, szptr, sxptr, syptr, rho, dt, dzi, dxi, dyi, 
         nz0, nzf, nx0, nxf, ny0, nyf, SZ, SX, SY, dimmz, dimmx);
#endif
};

#ifdef OPTIMIZED
template <const int HALO = 4, 
          const int BLOCK_DIM_X = 16,
          const int BLOCK_DIM_Y = 16>
__global__
__launch_bounds__(128, 16) 
void compute_component_vcell_BL_cuda_k ( float* __restrict__ vptr,
                                   const float* __restrict__ szptr,
                                   const float* __restrict__ sxptr,
                                   const float* __restrict__ syptr,
                                   const float* __restrict__ rho,
                                   const float           dt,
                                   const float           dzi,
                                   const float           dxi,
                                   const float           dyi,
                                   const int             nz0,
                                   const int             nzf,
                                   const int             nx0,
                                   const int             nxf,
                                   const int             ny0,
                                   const int             nyf,
                                   const int             SZ,
                                   const int             SX,
                                   const int             SY,
                                   const int             dimmz,
                                   const int             dimmx)
{
    for(int z = blockIdx.x * blockDim.x + threadIdx.x + nz0; 
            z < nzf; 
            z += gridDim.x * blockDim.x)
    {
        for(int x = blockIdx.y * blockDim.y + threadIdx.y + nx0; 
                x < nxf; 
                x += gridDim.y * blockDim.y)
        {
            __shared__ float sx_smem[BLOCK_DIM_Y][BLOCK_DIM_X+2*HALO];
            float sy_front1, sy_front2, sy_front3, sy_front4;
            float sy_back1, sy_back2, sy_back3, sy_back4;
            float sy_current;

            sy_back3   = syptr[IDX(z,x,ny0-HALO+0+SY,dimmz,dimmx)];
            sy_back2   = syptr[IDX(z,x,ny0-HALO+1+SY,dimmz,dimmx)];
            sy_back1   = syptr[IDX(z,x,ny0-HALO+2+SY,dimmz,dimmx)];
            sy_current = syptr[IDX(z,x,ny0-HALO+3+SY,dimmz,dimmx)];
            sy_front1  = syptr[IDX(z,x,ny0-HALO+4+SY,dimmz,dimmx)];
            sy_front2  = syptr[IDX(z,x,ny0-HALO+5+SY,dimmz,dimmx)];
            sy_front3  = syptr[IDX(z,x,ny0-HALO+6+SY,dimmz,dimmx)];
            sy_front4  = syptr[IDX(z,x,ny0-HALO+7+SY,dimmz,dimmx)];

            for(int y = ny0; 
                    y < nyf; 
                    y++)
            {
                /////////// register tiling-advance plane ////////////////
                sy_back4   = sy_back3;
                sy_back3   = sy_back2;
                sy_back2   = sy_back1;
                sy_back1   = sy_current;
                sy_current = sy_front1;
                sy_front1  = sy_front2;
                sy_front2  = sy_front3;
                sy_front3  = sy_front4;
                sy_front4  = syptr[IDX(z,x,y+SY+HALO,dimmz,dimmx)];
                //////////////////////////////////////////////////////////

                const int tx = threadIdx.x+HALO;
                const int ty = threadIdx.y;
                ///////////// intra-block communication///////////////////
                sx_smem[ty][tx] = sxptr[IDX(z,x+SX,y,dimmz,dimmx)];
                //if (threadIdx.y < HALO)
                //{
                //    sx_smem[threadIdx.y                 ][tx] = sxptr[IDX(z,x+SX-HALO,       y,dimmz,dimmx)];
                //    sx_smem[threadIdx.y+BLOCK_DIM_Y+HALO][tx] = sxptr[IDX(z,x+SX+BLOCK_DIM_X,y,dimmz,dimmx)];
                //}
                if (threadIdx.x < HALO)
                {
                    sx_smem[ty][threadIdx.x                 ] = sxptr[IDX(z,x+SX-HALO,       y,dimmz,dimmx)];
                    sx_smem[ty][threadIdx.x+BLOCK_DIM_X+HALO] = sxptr[IDX(z,x+SX+BLOCK_DIM_X,y,dimmz,dimmx)];
                }
                /////////////////////////////////////////////////////////
                
                ///////////// intra-warp communication /////////////////
                float sz_current = szptr[IDX(z+SZ,x,y,dimmz,dimmx)];
                float sz_front3 = __shfl_up(sz_current, 3);
                float sz_front2 = __shfl_up(sz_current, 2);
                float sz_front1 = __shfl_up(sz_current, 1);
                float sz_back1  = __shfl_down(sz_current, 1);
                float sz_back2  = __shfl_down(sz_current, 2);
                float sz_back3  = __shfl_down(sz_current, 3);
                float sz_back4  = __shfl_down(sz_current, 4);
                ////////////////////////////////////////////////////////
                float rho_current  = rho[IDX(z,x,y,dimmz,dimmx)];
                float rho_front1   = __shfl_up(rho_current, 1);
                ////////////////////////////////////////////////////////
                __syncthreads();

                const float lrho = (2.0f / (rho_current + rho_front1));

                const float stz = ((C0 * ( sz_current - sz_back1 ) +
                                    C1 * ( sz_front1  - sz_back2 ) +
                                    C2 * ( sz_front2  - sz_back3 ) +
                                    C3 * ( sz_front3  - sz_back4 )) * dzi );
                
                const float stx = ((C0 * ( sx_smem[ty][tx  ] - sx_smem[ty][tx-1] ) +
                                    C1 * ( sx_smem[ty][tx+1] - sx_smem[ty][tx-2] ) +
                                    C2 * ( sx_smem[ty][tx+2] - sx_smem[ty][tx-3] ) +
                                    C3 * ( sx_smem[ty][tx+3] - sx_smem[ty][tx-4] )) * dxi );

                const float sty = ((C0 * ( sy_current - sy_back1 ) +
                                    C1 * ( sy_front1  - sy_back2 ) +
                                    C2 * ( sy_front2  - sy_back3 ) +
                                    C3 * ( sy_front3  - sy_back4 )) * dyi );

                vptr[IDX(z,x,y,dimmz,dimmx)] += (stx  + sty  + stz) * dt * lrho;
            }
        }
    }
}
#else
__global__
void compute_component_vcell_BL_cuda_k ( float* __restrict__ vptr,
                                   const float* __restrict__ szptr,
                                   const float* __restrict__ sxptr,
                                   const float* __restrict__ syptr,
                                   const float* __restrict__ rho,
                                   const float           dt,
                                   const float           dzi,
                                   const float           dxi,
                                   const float           dyi,
                                   const int             nz0,
                                   const int             nzf,
                                   const int             nx0,
                                   const int             nxf,
                                   const int             ny0,
                                   const int             nyf,
                                   const int             SZ,
                                   const int             SX,
                                   const int             SY,
                                   const int             dimmz,
                                   const int             dimmx)
{
    for(int z = blockIdx.x * blockDim.x + threadIdx.x + nz0; 
            z < nzf; 
            z += gridDim.x * blockDim.x)
    {
        for(int x = blockIdx.y * blockDim.y + threadIdx.y + nx0; 
                x < nxf; 
                x += gridDim.y * blockDim.y)
        {
            for(int y = ny0; 
                    y < nyf; 
                    y++)
            {
                const float lrho = rho_BL(rho, z, x, y, dimmz, dimmx);
                
                const float stx  = stencil_X( SX, sxptr, dxi, z, x, y, dimmz, dimmx);
                const float sty  = stencil_Y( SY, syptr, dyi, z, x, y, dimmz, dimmx);
                const float stz  = stencil_Z( SZ, szptr, dzi, z, x, y, dimmz, dimmx);
                
                vptr[IDX(z,x,y,dimmz,dimmx)] += (stx  + sty  + stz) * dt * lrho;
            }
        }
    }
}
#endif

extern "C"
void compute_component_vcell_BL_cuda ( float* vptr,
                                 const float* szptr,
                                 const float* sxptr,
                                 const float* syptr,
                                 const float* rho,
                                 const float  dt,
                                 const float  dzi,
                                 const float  dxi,
                                 const float  dyi,
                                 const int    nz0,
                                 const int    nzf,
                                 const int    nx0,
                                 const int    nxf,
                                 const int    ny0,
                                 const int    nyf,
                                 const int    SZ,
                                 const int    SX,
                                 const int    SY,
                                 const int    dimmz,
                                 const int    dimmx,
                                 void*        stream)
{
    const int block_dim_x = 32;
    const int block_dim_y = 4;


    dim3 grid_dim( ((nzf-nz0) + block_dim_x-1)/block_dim_x,
                   ((nxf-nx0) + block_dim_y-1)/block_dim_y,
                      1 );
    dim3 block_dim(block_dim_x, block_dim_y, 1);

    cudaStream_t s = (cudaStream_t) stream;

#ifdef OPTIMIZED
    compute_component_vcell_BL_cuda_k<4,block_dim_x,block_dim_y><<<grid_dim, block_dim, 0, s>>>
        (vptr, szptr, sxptr, syptr, rho, dt, dzi, dxi, dyi, 
         nz0, nzf, nx0, nxf, ny0, nyf, SZ, SX, SY, dimmz, dimmx);
#else
    compute_component_vcell_BL_cuda_k<<<grid_dim, block_dim, 0, s>>>
        (vptr, szptr, sxptr, syptr, rho, dt, dzi, dxi, dyi, 
         nz0, nzf, nx0, nxf, ny0, nyf, SZ, SX, SY, dimmz, dimmx);
#endif
};





/* ------------------------------------------------------------------------------ */
/*                                                                                */
/*                               CALCULO DE TENSIONES                             */
/*                                                                                */
/* ------------------------------------------------------------------------------ */

__device__ inline
void stress_update(float* __restrict__ sptr,
                   const float     c1,
                   const float     c2,
                   const float     c3,
                   const float     c4,
                   const float     c5,
                   const float     c6,
                   const int  z,
                   const int  x,
                   const int  y,
                   const float     dt,
                   const float     u_x,
                   const float     u_y,
                   const float     u_z,
                   const float     v_x,
                   const float     v_y,
                   const float     v_z,
                   const float     w_x,
                   const float     w_y,
                   const float     w_z,
                   const int  dimmz,
                   const int  dimmx)
{
    float accum  = dt * c1 * u_x;
         accum += dt * c2 * v_y;
         accum += dt * c3 * w_z;
         accum += dt * c4 * (w_y + v_z);
         accum += dt * c5 * (w_x + u_z);
         accum += dt * c6 * (v_x + u_y);
    sptr[IDX(z,x,y,dimmz,dimmx)] += accum;
};

__device__ inline
float cell_coeff_BR (const float* __restrict__ ptr, 
                     const int z,
                     const int x,
                     const int y,
                     const int dimmz,
                     const int dimmx)
{
    return ( 1.0f / ( 2.5f  *(ptr[IDX(z  , x  ,y,dimmz,dimmx)] +
                              ptr[IDX(z  , x+1,y,dimmz,dimmx)] +
                              ptr[IDX(z+1, x  ,y,dimmz,dimmx)] +
                              ptr[IDX(z+1, x+1,y,dimmz,dimmx)])) );
};

__device__ inline
float cell_coeff_TL (const float* __restrict__ ptr, 
                     const int z, 
                     const int x, 
                     const int y, 
                     const int dimmz, 
                     const int dimmx)
{
    return ( 1.0f / (ptr[IDX(z,x,y,dimmz,dimmx)]));
};

__device__ inline
float cell_coeff_BL (const float* __restrict__ ptr, 
                     const int z, 
                     const int x, 
                     const int y, 
                     const int dimmz, 
                     const int dimmx)
{
    return ( 1.0f / ( 2.5f *(ptr[IDX(z  ,x,y  ,dimmz,dimmx)] +
                             ptr[IDX(z  ,x,y+1,dimmz,dimmx)] +
                             ptr[IDX(z+1,x,y  ,dimmz,dimmx)] +
                             ptr[IDX(z+1,x,y+1,dimmz,dimmx)])) );
};

__device__ inline
float cell_coeff_TR (const float* __restrict__ ptr, 
                     const int z, 
                     const int x, 
                     const int y, 
                     const int dimmz, 
                     const int dimmx)
{
    return ( 1.0f / ( 2.5f *(ptr[IDX(z  , x  , y  ,dimmz,dimmx)] +
                             ptr[IDX(z  , x+1, y  ,dimmz,dimmx)] +
                             ptr[IDX(z  , x  , y+1,dimmz,dimmx)] +
                             ptr[IDX(z  , x+1, y+1,dimmz,dimmx)])));
};

__device__ inline
float cell_coeff_ARTM_BR(const float* __restrict__ ptr, 
                         const int z, 
                         const int x, 
                         const int y, 
                         const int dimmz, 
                         const int dimmx)
{
    return ((1.0f / ptr[IDX(z  ,x  ,y,dimmz,dimmx )]  +
             1.0f / ptr[IDX(z  ,x+1,y,dimmz,dimmx )]  +
             1.0f / ptr[IDX(z+1,x  ,y,dimmz,dimmx )]  +
             1.0f / ptr[IDX(z+1,x+1,y,dimmz,dimmx )]) * 0.25f);
};

__device__ inline
float cell_coeff_ARTM_TL( const float* __restrict__ ptr, 
                         const int z, 
                         const int x, 
                         const int y, 
                         const int dimmz, 
                         const int dimmx)
{
    return (1.0f / ptr[IDX(z,x,y,dimmz,dimmx)]);
};

__device__ inline
float cell_coeff_ARTM_BL(const float* __restrict__ ptr, 
                         const int z, 
                         const int x, 
                         const int y, 
                         const int dimmz, 
                         const int dimmx)
{
    return ((1.0f / ptr[IDX(z  ,x,y  ,dimmz,dimmx)]  +
             1.0f / ptr[IDX(z  ,x,y+1,dimmz,dimmx)]  +
             1.0f / ptr[IDX(z+1,x,y  ,dimmz,dimmx)]  +
             1.0f / ptr[IDX(z+1,x,y+1,dimmz,dimmx)]) * 0.25f);
};

__device__ inline
float cell_coeff_ARTM_TR(const float* __restrict__ ptr, 
                         const int z, 
                         const int x, 
                         const int y, 
                         const int dimmz, 
                         const int dimmx)
{
    return ((1.0f / ptr[IDX(z,x  ,y  ,dimmz,dimmx)]  +
             1.0f / ptr[IDX(z,x+1,y  ,dimmz,dimmx)]  +
             1.0f / ptr[IDX(z,x  ,y+1,dimmz,dimmx)]  +
             1.0f / ptr[IDX(z,x+1,y+1,dimmz,dimmx)]) * 0.25f);
};



__global__
__launch_bounds__(128, 8) 
void compute_component_scell_TR_cuda_k ( float* __restrict__ sxxptr,
                                         float* __restrict__ syyptr,
                                         float* __restrict__ szzptr,
                                         float* __restrict__ syzptr,
                                         float* __restrict__ sxzptr,
                                         float* __restrict__ sxyptr,
                                   const float* __restrict__ vxu,
                                   const float* __restrict__ vxv,
                                   const float* __restrict__ vxw,
                                   const float* __restrict__ vyu,
                                   const float* __restrict__ vyv,
                                   const float* __restrict__ vyw,
                                   const float* __restrict__ vzu,
                                   const float* __restrict__ vzv,
                                   const float* __restrict__ vzw,
                                   const float* __restrict__ cc11,
                                   const float* __restrict__ cc12,
                                   const float* __restrict__ cc13,
                                   const float* __restrict__ cc14,
                                   const float* __restrict__ cc15,
                                   const float* __restrict__ cc16,
                                   const float* __restrict__ cc22,
                                   const float* __restrict__ cc23,
                                   const float* __restrict__ cc24,
                                   const float* __restrict__ cc25,
                                   const float* __restrict__ cc26,
                                   const float* __restrict__ cc33,
                                   const float* __restrict__ cc34,
                                   const float* __restrict__ cc35,
                                   const float* __restrict__ cc36,
                                   const float* __restrict__ cc44,
                                   const float* __restrict__ cc45,
                                   const float* __restrict__ cc46,
                                   const float* __restrict__ cc55,
                                   const float* __restrict__ cc56,
                                   const float* __restrict__ cc66,
                                   const float  dt,
                                   const float  dzi,
                                   const float  dxi,
                                   const float  dyi,
                                   const int    nz0,
                                   const int    nzf,
                                   const int    nx0,
                                   const int    nxf,
                                   const int    ny0,
                                   const int    nyf,
                                   const int    SZ,
                                   const int    SX,
                                   const int    SY,
                                   const int    dimmz,
                                   const int    dimmx)
{
    for(int z = blockIdx.x * blockDim.x + threadIdx.x + nz0; 
            z < nzf; 
            z += gridDim.x * blockDim.x)
    {
        for(int x = blockIdx.y * blockDim.y + threadIdx.y + nx0; 
                x < nxf; 
                x += gridDim.y * blockDim.y)
        {
            for(int y = ny0; 
                    y < nyf; 
                    y++)
            {
                const float c11 = cell_coeff_TR      (cc11, z, x, y, dimmz, dimmx);
                const float c12 = cell_coeff_TR      (cc12, z, x, y, dimmz, dimmx);
                const float c13 = cell_coeff_TR      (cc13, z, x, y, dimmz, dimmx);
                const float c14 = cell_coeff_ARTM_TR (cc14, z, x, y, dimmz, dimmx);
                const float c15 = cell_coeff_ARTM_TR (cc15, z, x, y, dimmz, dimmx);
                const float c16 = cell_coeff_ARTM_TR (cc16, z, x, y, dimmz, dimmx);
                const float c22 = cell_coeff_TR      (cc22, z, x, y, dimmz, dimmx);
                const float c23 = cell_coeff_TR      (cc23, z, x, y, dimmz, dimmx);
                const float c24 = cell_coeff_ARTM_TR (cc24, z, x, y, dimmz, dimmx);
                const float c25 = cell_coeff_ARTM_TR (cc25, z, x, y, dimmz, dimmx);
                const float c26 = cell_coeff_ARTM_TR (cc26, z, x, y, dimmz, dimmx);
                const float c33 = cell_coeff_TR      (cc33, z, x, y, dimmz, dimmx);
                const float c34 = cell_coeff_ARTM_TR (cc34, z, x, y, dimmz, dimmx);
                const float c35 = cell_coeff_ARTM_TR (cc35, z, x, y, dimmz, dimmx);
                const float c36 = cell_coeff_ARTM_TR (cc36, z, x, y, dimmz, dimmx);
                const float c44 = cell_coeff_TR      (cc44, z, x, y, dimmz, dimmx);
                const float c45 = cell_coeff_ARTM_TR (cc45, z, x, y, dimmz, dimmx);
                const float c46 = cell_coeff_ARTM_TR (cc46, z, x, y, dimmz, dimmx);
                const float c55 = cell_coeff_TR      (cc55, z, x, y, dimmz, dimmx);
                const float c56 = cell_coeff_ARTM_TR (cc56, z, x, y, dimmz, dimmx);
                const float c66 = cell_coeff_TR      (cc66, z, x, y, dimmz, dimmx);
                
                const float u_x = stencil_X (SX, vxu, dxi, z, x, y, dimmz, dimmx);
                const float v_x = stencil_X (SX, vxv, dxi, z, x, y, dimmz, dimmx);
                const float w_x = stencil_X (SX, vxw, dxi, z, x, y, dimmz, dimmx);
                
                const float u_y = stencil_Y (SY, vyu, dyi, z, x, y, dimmz, dimmx);
                const float v_y = stencil_Y (SY, vyv, dyi, z, x, y, dimmz, dimmx);
                const float w_y = stencil_Y (SY, vyw, dyi, z, x, y, dimmz, dimmx);
                
                const float u_z = stencil_Z (SZ, vzu, dzi, z, x, y, dimmz, dimmx);
                const float v_z = stencil_Z (SZ, vzv, dzi, z, x, y, dimmz, dimmx);
                const float w_z = stencil_Z (SZ, vzw, dzi, z, x, y, dimmz, dimmx);
                
                stress_update (sxxptr,c11,c12,c13,c14,c15,c16,z,x,y,dt,u_x,u_y,u_z,v_x,v_y,v_z,w_x,w_y,w_z,dimmz,dimmx );
                stress_update (syyptr,c12,c22,c23,c24,c25,c26,z,x,y,dt,u_x,u_y,u_z,v_x,v_y,v_z,w_x,w_y,w_z,dimmz,dimmx );
                stress_update (szzptr,c13,c23,c33,c34,c35,c36,z,x,y,dt,u_x,u_y,u_z,v_x,v_y,v_z,w_x,w_y,w_z,dimmz,dimmx );
                stress_update (syzptr,c14,c24,c34,c44,c45,c46,z,x,y,dt,u_x,u_y,u_z,v_x,v_y,v_z,w_x,w_y,w_z,dimmz,dimmx );
                stress_update (sxzptr,c15,c25,c35,c45,c55,c56,z,x,y,dt,u_x,u_y,u_z,v_x,v_y,v_z,w_x,w_y,w_z,dimmz,dimmx );
                stress_update (sxyptr,c16,c26,c36,c46,c56,c66,z,x,y,dt,u_x,u_y,u_z,v_x,v_y,v_z,w_x,w_y,w_z,dimmz,dimmx );
            }
        }
    }
}

extern "C"
void compute_component_scell_TR_cuda ( float* sxxptr,
                                       float* syyptr,
                                       float* szzptr,
                                       float* syzptr,
                                       float* sxzptr,
                                       float* sxyptr,
                                 const float* vxu,
                                 const float* vxv,
                                 const float* vxw,
                                 const float* vyu,
                                 const float* vyv,
                                 const float* vyw,
                                 const float* vzu,
                                 const float* vzv,
                                 const float* vzw,
                                 const float* cc11,
                                 const float* cc12,
                                 const float* cc13,
                                 const float* cc14,
                                 const float* cc15,
                                 const float* cc16,
                                 const float* cc22,
                                 const float* cc23,
                                 const float* cc24,
                                 const float* cc25,
                                 const float* cc26,
                                 const float* cc33,
                                 const float* cc34,
                                 const float* cc35,
                                 const float* cc36,
                                 const float* cc44,
                                 const float* cc45,
                                 const float* cc46,
                                 const float* cc55,
                                 const float* cc56,
                                 const float* cc66,
                                 const float  dt,
                                 const float  dzi,
                                 const float  dxi,
                                 const float  dyi,
                                 const int    nz0,
                                 const int    nzf,
                                 const int    nx0,
                                 const int    nxf,
                                 const int    ny0,
                                 const int    nyf,
                                 const int    SZ,
                                 const int    SX,
                                 const int    SY,
                                 const int    dimmz,
                                 const int    dimmx,
                                 void*        stream)
{
    const int block_dim_x = 32;
    const int block_dim_y = 4;


    dim3 grid_dim( ((nzf-nz0) + block_dim_x-1)/block_dim_x,
                   ((nxf-nx0) + block_dim_y-1)/block_dim_y,
                      1 );
    dim3 block_dim(block_dim_x, block_dim_y, 1);

    cudaStream_t s = (cudaStream_t) stream;

    compute_component_scell_TR_cuda_k<<<grid_dim, block_dim, 0, s>>>
        (sxxptr, syyptr, szzptr, syzptr, sxzptr, sxyptr,
         vxu, vxv, vxw, vyu, vyv, vyw, vzu, vzv, vzw,
         cc11, cc12, cc13, cc14, cc15, cc16,
         cc22, cc23, cc24, cc25, cc26,
         cc33, cc34, cc35, cc36,
         cc44, cc45, cc45,
         cc55, cc56,
         cc66,
         dt, dzi, dxi, dyi, nz0, nzf, nx0, nxf, ny0, nyf, SZ, SX, SY, dimmz, dimmx);
};


__global__
__launch_bounds__(128, 8) 
void compute_component_scell_TL_cuda_k ( float* __restrict__ sxxptr,
                                         float* __restrict__ syyptr,
                                         float* __restrict__ szzptr,
                                         float* __restrict__ syzptr,
                                         float* __restrict__ sxzptr,
                                         float* __restrict__ sxyptr,
                                   const float* __restrict__ vxu,
                                   const float* __restrict__ vxv,
                                   const float* __restrict__ vxw,
                                   const float* __restrict__ vyu,
                                   const float* __restrict__ vyv,
                                   const float* __restrict__ vyw,
                                   const float* __restrict__ vzu,
                                   const float* __restrict__ vzv,
                                   const float* __restrict__ vzw,
                                   const float* __restrict__ cc11,
                                   const float* __restrict__ cc12,
                                   const float* __restrict__ cc13,
                                   const float* __restrict__ cc14,
                                   const float* __restrict__ cc15,
                                   const float* __restrict__ cc16,
                                   const float* __restrict__ cc22,
                                   const float* __restrict__ cc23,
                                   const float* __restrict__ cc24,
                                   const float* __restrict__ cc25,
                                   const float* __restrict__ cc26,
                                   const float* __restrict__ cc33,
                                   const float* __restrict__ cc34,
                                   const float* __restrict__ cc35,
                                   const float* __restrict__ cc36,
                                   const float* __restrict__ cc44,
                                   const float* __restrict__ cc45,
                                   const float* __restrict__ cc46,
                                   const float* __restrict__ cc55,
                                   const float* __restrict__ cc56,
                                   const float* __restrict__ cc66,
                                   const float  dt,
                                   const float  dzi,
                                   const float  dxi,
                                   const float  dyi,
                                   const int    nz0,
                                   const int    nzf,
                                   const int    nx0,
                                   const int    nxf,
                                   const int    ny0,
                                   const int    nyf,
                                   const int    SZ,
                                   const int    SX,
                                   const int    SY,
                                   const int    dimmz,
                                   const int    dimmx)
{
    for(int z = blockIdx.x * blockDim.x + threadIdx.x + nz0; 
            z < nzf; 
            z += gridDim.x * blockDim.x)
    {
        for(int x = blockIdx.y * blockDim.y + threadIdx.y + nx0; 
                x < nxf; 
                x += gridDim.y * blockDim.y)
        {
            for(int y = ny0; 
                    y < nyf; 
                    y++)
            {
                const float c11 = cell_coeff_TL      (cc11, z, x, y, dimmz, dimmx);
                const float c12 = cell_coeff_TL      (cc12, z, x, y, dimmz, dimmx);
                const float c13 = cell_coeff_TL      (cc13, z, x, y, dimmz, dimmx);
                const float c14 = cell_coeff_ARTM_TL (cc14, z, x, y, dimmz, dimmx);
                const float c15 = cell_coeff_ARTM_TL (cc15, z, x, y, dimmz, dimmx);
                const float c16 = cell_coeff_ARTM_TL (cc16, z, x, y, dimmz, dimmx);
                const float c22 = cell_coeff_TL      (cc22, z, x, y, dimmz, dimmx);
                const float c23 = cell_coeff_TL      (cc23, z, x, y, dimmz, dimmx);
                const float c24 = cell_coeff_ARTM_TL (cc24, z, x, y, dimmz, dimmx);
                const float c25 = cell_coeff_ARTM_TL (cc25, z, x, y, dimmz, dimmx);
                const float c26 = cell_coeff_ARTM_TL (cc26, z, x, y, dimmz, dimmx);
                const float c33 = cell_coeff_TL      (cc33, z, x, y, dimmz, dimmx);
                const float c34 = cell_coeff_ARTM_TL (cc34, z, x, y, dimmz, dimmx);
                const float c35 = cell_coeff_ARTM_TL (cc35, z, x, y, dimmz, dimmx);
                const float c36 = cell_coeff_ARTM_TL (cc36, z, x, y, dimmz, dimmx);
                const float c44 = cell_coeff_TL      (cc44, z, x, y, dimmz, dimmx);
                const float c45 = cell_coeff_ARTM_TL (cc45, z, x, y, dimmz, dimmx);
                const float c46 = cell_coeff_ARTM_TL (cc46, z, x, y, dimmz, dimmx);
                const float c55 = cell_coeff_TL      (cc55, z, x, y, dimmz, dimmx);
                const float c56 = cell_coeff_ARTM_TL (cc56, z, x, y, dimmz, dimmx);
                const float c66 = cell_coeff_TL      (cc66, z, x, y, dimmz, dimmx);
                
                const float u_x = stencil_X (SX, vxu, dxi, z, x, y, dimmz, dimmx);
                const float v_x = stencil_X (SX, vxv, dxi, z, x, y, dimmz, dimmx);
                const float w_x = stencil_X (SX, vxw, dxi, z, x, y, dimmz, dimmx);
                
                const float u_y = stencil_Y (SY, vyu, dyi, z, x, y, dimmz, dimmx);
                const float v_y = stencil_Y (SY, vyv, dyi, z, x, y, dimmz, dimmx);
                const float w_y = stencil_Y (SY, vyw, dyi, z, x, y, dimmz, dimmx);
                
                const float u_z = stencil_Z (SZ, vzu, dzi, z, x, y, dimmz, dimmx);
                const float v_z = stencil_Z (SZ, vzv, dzi, z, x, y, dimmz, dimmx);
                const float w_z = stencil_Z (SZ, vzw, dzi, z, x, y, dimmz, dimmx);
                
                stress_update (sxxptr,c11,c12,c13,c14,c15,c16,z,x,y,dt,u_x,u_y,u_z,v_x,v_y,v_z,w_x,w_y,w_z,dimmz,dimmx );
                stress_update (syyptr,c12,c22,c23,c24,c25,c26,z,x,y,dt,u_x,u_y,u_z,v_x,v_y,v_z,w_x,w_y,w_z,dimmz,dimmx );
                stress_update (szzptr,c13,c23,c33,c34,c35,c36,z,x,y,dt,u_x,u_y,u_z,v_x,v_y,v_z,w_x,w_y,w_z,dimmz,dimmx );
                stress_update (syzptr,c14,c24,c34,c44,c45,c46,z,x,y,dt,u_x,u_y,u_z,v_x,v_y,v_z,w_x,w_y,w_z,dimmz,dimmx );
                stress_update (sxzptr,c15,c25,c35,c45,c55,c56,z,x,y,dt,u_x,u_y,u_z,v_x,v_y,v_z,w_x,w_y,w_z,dimmz,dimmx );
                stress_update (sxyptr,c16,c26,c36,c46,c56,c66,z,x,y,dt,u_x,u_y,u_z,v_x,v_y,v_z,w_x,w_y,w_z,dimmz,dimmx );
            }
        }
    }
}

extern "C"
void compute_component_scell_TL_cuda ( float* sxxptr,
                                       float* syyptr,
                                       float* szzptr,
                                       float* syzptr,
                                       float* sxzptr,
                                       float* sxyptr,
                                 const float* vxu,
                                 const float* vxv,
                                 const float* vxw,
                                 const float* vyu,
                                 const float* vyv,
                                 const float* vyw,
                                 const float* vzu,
                                 const float* vzv,
                                 const float* vzw,
                                 const float* cc11,
                                 const float* cc12,
                                 const float* cc13,
                                 const float* cc14,
                                 const float* cc15,
                                 const float* cc16,
                                 const float* cc22,
                                 const float* cc23,
                                 const float* cc24,
                                 const float* cc25,
                                 const float* cc26,
                                 const float* cc33,
                                 const float* cc34,
                                 const float* cc35,
                                 const float* cc36,
                                 const float* cc44,
                                 const float* cc45,
                                 const float* cc46,
                                 const float* cc55,
                                 const float* cc56,
                                 const float* cc66,
                                 const float  dt,
                                 const float  dzi,
                                 const float  dxi,
                                 const float  dyi,
                                 const int    nz0,
                                 const int    nzf,
                                 const int    nx0,
                                 const int    nxf,
                                 const int    ny0,
                                 const int    nyf,
                                 const int    SZ,
                                 const int    SX,
                                 const int    SY,
                                 const int    dimmz,
                                 const int    dimmx,
                                 void*        stream)
{
    const int block_dim_x = 32;
    const int block_dim_y = 4;


    dim3 grid_dim( ((nzf-nz0) + block_dim_x-1)/block_dim_x,
                   ((nxf-nx0) + block_dim_y-1)/block_dim_y,
                      1 );
    dim3 block_dim(block_dim_x, block_dim_y, 1);

    cudaStream_t s = (cudaStream_t) stream;

    compute_component_scell_TL_cuda_k<<<grid_dim, block_dim, 0, s>>>
        (sxxptr, syyptr, szzptr, syzptr, sxzptr, sxyptr,
         vxu, vxv, vxw, vyu, vyv, vyw, vzu, vzv, vzw,
         cc11, cc12, cc13, cc14, cc15, cc16,
         cc22, cc23, cc24, cc25, cc26,
         cc33, cc34, cc35, cc36,
         cc44, cc45, cc45,
         cc55, cc56,
         cc66,
         dt, dzi, dxi, dyi, nz0, nzf, nx0, nxf, ny0, nyf, SZ, SX, SY, dimmz, dimmx);
};


__global__
__launch_bounds__(128, 8) 
void compute_component_scell_BR_cuda_k ( float* __restrict__ sxxptr,
                                         float* __restrict__ syyptr,
                                         float* __restrict__ szzptr,
                                         float* __restrict__ syzptr,
                                         float* __restrict__ sxzptr,
                                         float* __restrict__ sxyptr,
                                   const float* __restrict__ vxu,
                                   const float* __restrict__ vxv,
                                   const float* __restrict__ vxw,
                                   const float* __restrict__ vyu,
                                   const float* __restrict__ vyv,
                                   const float* __restrict__ vyw,
                                   const float* __restrict__ vzu,
                                   const float* __restrict__ vzv,
                                   const float* __restrict__ vzw,
                                   const float* __restrict__ cc11,
                                   const float* __restrict__ cc12,
                                   const float* __restrict__ cc13,
                                   const float* __restrict__ cc14,
                                   const float* __restrict__ cc15,
                                   const float* __restrict__ cc16,
                                   const float* __restrict__ cc22,
                                   const float* __restrict__ cc23,
                                   const float* __restrict__ cc24,
                                   const float* __restrict__ cc25,
                                   const float* __restrict__ cc26,
                                   const float* __restrict__ cc33,
                                   const float* __restrict__ cc34,
                                   const float* __restrict__ cc35,
                                   const float* __restrict__ cc36,
                                   const float* __restrict__ cc44,
                                   const float* __restrict__ cc45,
                                   const float* __restrict__ cc46,
                                   const float* __restrict__ cc55,
                                   const float* __restrict__ cc56,
                                   const float* __restrict__ cc66,
                                   const float  dt,
                                   const float  dzi,
                                   const float  dxi,
                                   const float  dyi,
                                   const int    nz0,
                                   const int    nzf,
                                   const int    nx0,
                                   const int    nxf,
                                   const int    ny0,
                                   const int    nyf,
                                   const int    SZ,
                                   const int    SX,
                                   const int    SY,
                                   const int    dimmz,
                                   const int    dimmx)
{
    for(int z = blockIdx.x * blockDim.x + threadIdx.x + nz0; 
            z < nzf; 
            z += gridDim.x * blockDim.x)
    {
        for(int x = blockIdx.y * blockDim.y + threadIdx.y + nx0; 
                x < nxf; 
                x += gridDim.y * blockDim.y)
        {
            for(int y = ny0; 
                    y < nyf; 
                    y++)
            {
                const float c11 = cell_coeff_BR      (cc11, z, x, y, dimmz, dimmx);
                const float c12 = cell_coeff_BR      (cc12, z, x, y, dimmz, dimmx);
                const float c13 = cell_coeff_BR      (cc13, z, x, y, dimmz, dimmx);
                const float c14 = cell_coeff_ARTM_BR (cc14, z, x, y, dimmz, dimmx);
                const float c15 = cell_coeff_ARTM_BR (cc15, z, x, y, dimmz, dimmx);
                const float c16 = cell_coeff_ARTM_BR (cc16, z, x, y, dimmz, dimmx);
                const float c22 = cell_coeff_BR      (cc22, z, x, y, dimmz, dimmx);
                const float c23 = cell_coeff_BR      (cc23, z, x, y, dimmz, dimmx);
                const float c24 = cell_coeff_ARTM_BR (cc24, z, x, y, dimmz, dimmx);
                const float c25 = cell_coeff_ARTM_BR (cc25, z, x, y, dimmz, dimmx);
                const float c26 = cell_coeff_ARTM_BR (cc26, z, x, y, dimmz, dimmx);
                const float c33 = cell_coeff_BR      (cc33, z, x, y, dimmz, dimmx);
                const float c34 = cell_coeff_ARTM_BR (cc34, z, x, y, dimmz, dimmx);
                const float c35 = cell_coeff_ARTM_BR (cc35, z, x, y, dimmz, dimmx);
                const float c36 = cell_coeff_ARTM_BR (cc36, z, x, y, dimmz, dimmx);
                const float c44 = cell_coeff_BR      (cc44, z, x, y, dimmz, dimmx);
                const float c45 = cell_coeff_ARTM_BR (cc45, z, x, y, dimmz, dimmx);
                const float c46 = cell_coeff_ARTM_BR (cc46, z, x, y, dimmz, dimmx);
                const float c55 = cell_coeff_BR      (cc55, z, x, y, dimmz, dimmx);
                const float c56 = cell_coeff_ARTM_BR (cc56, z, x, y, dimmz, dimmx);
                const float c66 = cell_coeff_BR      (cc66, z, x, y, dimmz, dimmx);
                
                const float u_x = stencil_X (SX, vxu, dxi, z, x, y, dimmz, dimmx);
                const float v_x = stencil_X (SX, vxv, dxi, z, x, y, dimmz, dimmx);
                const float w_x = stencil_X (SX, vxw, dxi, z, x, y, dimmz, dimmx);
                
                const float u_y = stencil_Y (SY, vyu, dyi, z, x, y, dimmz, dimmx);
                const float v_y = stencil_Y (SY, vyv, dyi, z, x, y, dimmz, dimmx);
                const float w_y = stencil_Y (SY, vyw, dyi, z, x, y, dimmz, dimmx);
                
                const float u_z = stencil_Z (SZ, vzu, dzi, z, x, y, dimmz, dimmx);
                const float v_z = stencil_Z (SZ, vzv, dzi, z, x, y, dimmz, dimmx);
                const float w_z = stencil_Z (SZ, vzw, dzi, z, x, y, dimmz, dimmx);
                
                stress_update (sxxptr,c11,c12,c13,c14,c15,c16,z,x,y,dt,u_x,u_y,u_z,v_x,v_y,v_z,w_x,w_y,w_z,dimmz,dimmx );
                stress_update (syyptr,c12,c22,c23,c24,c25,c26,z,x,y,dt,u_x,u_y,u_z,v_x,v_y,v_z,w_x,w_y,w_z,dimmz,dimmx );
                stress_update (szzptr,c13,c23,c33,c34,c35,c36,z,x,y,dt,u_x,u_y,u_z,v_x,v_y,v_z,w_x,w_y,w_z,dimmz,dimmx );
                stress_update (syzptr,c14,c24,c34,c44,c45,c46,z,x,y,dt,u_x,u_y,u_z,v_x,v_y,v_z,w_x,w_y,w_z,dimmz,dimmx );
                stress_update (sxzptr,c15,c25,c35,c45,c55,c56,z,x,y,dt,u_x,u_y,u_z,v_x,v_y,v_z,w_x,w_y,w_z,dimmz,dimmx );
                stress_update (sxyptr,c16,c26,c36,c46,c56,c66,z,x,y,dt,u_x,u_y,u_z,v_x,v_y,v_z,w_x,w_y,w_z,dimmz,dimmx );
            }
        }
    }
}

extern "C"
void compute_component_scell_BR_cuda ( float* sxxptr,
                                       float* syyptr,
                                       float* szzptr,
                                       float* syzptr,
                                       float* sxzptr,
                                       float* sxyptr,
                                 const float* vxu,
                                 const float* vxv,
                                 const float* vxw,
                                 const float* vyu,
                                 const float* vyv,
                                 const float* vyw,
                                 const float* vzu,
                                 const float* vzv,
                                 const float* vzw,
                                 const float* cc11,
                                 const float* cc12,
                                 const float* cc13,
                                 const float* cc14,
                                 const float* cc15,
                                 const float* cc16,
                                 const float* cc22,
                                 const float* cc23,
                                 const float* cc24,
                                 const float* cc25,
                                 const float* cc26,
                                 const float* cc33,
                                 const float* cc34,
                                 const float* cc35,
                                 const float* cc36,
                                 const float* cc44,
                                 const float* cc45,
                                 const float* cc46,
                                 const float* cc55,
                                 const float* cc56,
                                 const float* cc66,
                                 const float  dt,
                                 const float  dzi,
                                 const float  dxi,
                                 const float  dyi,
                                 const int    nz0,
                                 const int    nzf,
                                 const int    nx0,
                                 const int    nxf,
                                 const int    ny0,
                                 const int    nyf,
                                 const int    SZ,
                                 const int    SX,
                                 const int    SY,
                                 const int    dimmz,
                                 const int    dimmx,
                                 void*        stream)
{
    const int block_dim_x = 32;
    const int block_dim_y = 4;


    dim3 grid_dim( ((nzf-nz0) + block_dim_x-1)/block_dim_x,
                   ((nxf-nx0) + block_dim_y-1)/block_dim_y,
                      1 );
    dim3 block_dim(block_dim_x, block_dim_y, 1);

    cudaStream_t s = (cudaStream_t) stream;

    compute_component_scell_BR_cuda_k<<<grid_dim, block_dim, 0, s>>>
        (sxxptr, syyptr, szzptr, syzptr, sxzptr, sxyptr,
         vxu, vxv, vxw, vyu, vyv, vyw, vzu, vzv, vzw,
         cc11, cc12, cc13, cc14, cc15, cc16,
         cc22, cc23, cc24, cc25, cc26,
         cc33, cc34, cc35, cc36,
         cc44, cc45, cc45,
         cc55, cc56,
         cc66,
         dt, dzi, dxi, dyi, nz0, nzf, nx0, nxf, ny0, nyf, SZ, SX, SY, dimmz, dimmx);
};


__global__
__launch_bounds__(128, 8) 
void compute_component_scell_BL_cuda_k ( float* __restrict__ sxxptr,
                                         float* __restrict__ syyptr,
                                         float* __restrict__ szzptr,
                                         float* __restrict__ syzptr,
                                         float* __restrict__ sxzptr,
                                         float* __restrict__ sxyptr,
                                   const float* __restrict__ vxu,
                                   const float* __restrict__ vxv,
                                   const float* __restrict__ vxw,
                                   const float* __restrict__ vyu,
                                   const float* __restrict__ vyv,
                                   const float* __restrict__ vyw,
                                   const float* __restrict__ vzu,
                                   const float* __restrict__ vzv,
                                   const float* __restrict__ vzw,
                                   const float* __restrict__ cc11,
                                   const float* __restrict__ cc12,
                                   const float* __restrict__ cc13,
                                   const float* __restrict__ cc14,
                                   const float* __restrict__ cc15,
                                   const float* __restrict__ cc16,
                                   const float* __restrict__ cc22,
                                   const float* __restrict__ cc23,
                                   const float* __restrict__ cc24,
                                   const float* __restrict__ cc25,
                                   const float* __restrict__ cc26,
                                   const float* __restrict__ cc33,
                                   const float* __restrict__ cc34,
                                   const float* __restrict__ cc35,
                                   const float* __restrict__ cc36,
                                   const float* __restrict__ cc44,
                                   const float* __restrict__ cc45,
                                   const float* __restrict__ cc46,
                                   const float* __restrict__ cc55,
                                   const float* __restrict__ cc56,
                                   const float* __restrict__ cc66,
                                   const float  dt,
                                   const float  dzi,
                                   const float  dxi,
                                   const float  dyi,
                                   const int    nz0,
                                   const int    nzf,
                                   const int    nx0,
                                   const int    nxf,
                                   const int    ny0,
                                   const int    nyf,
                                   const int    SZ,
                                   const int    SX,
                                   const int    SY,
                                   const int    dimmz,
                                   const int    dimmx)
{
    for(int z = blockIdx.x * blockDim.x + threadIdx.x + nz0; 
            z < nzf; 
            z += gridDim.x * blockDim.x)
    {
        for(int x = blockIdx.y * blockDim.y + threadIdx.y + nx0; 
                x < nxf; 
                x += gridDim.y * blockDim.y)
        {
            for(int y = ny0; 
                    y < nyf; 
                    y++)
            {
                const float c11 = cell_coeff_BL      (cc11, z, x, y, dimmz, dimmx);
                const float c12 = cell_coeff_BL      (cc12, z, x, y, dimmz, dimmx);
                const float c13 = cell_coeff_BL      (cc13, z, x, y, dimmz, dimmx);
                const float c14 = cell_coeff_ARTM_BL (cc14, z, x, y, dimmz, dimmx);
                const float c15 = cell_coeff_ARTM_BL (cc15, z, x, y, dimmz, dimmx);
                const float c16 = cell_coeff_ARTM_BL (cc16, z, x, y, dimmz, dimmx);
                const float c22 = cell_coeff_BL      (cc22, z, x, y, dimmz, dimmx);
                const float c23 = cell_coeff_BL      (cc23, z, x, y, dimmz, dimmx);
                const float c24 = cell_coeff_ARTM_BL (cc24, z, x, y, dimmz, dimmx);
                const float c25 = cell_coeff_ARTM_BL (cc25, z, x, y, dimmz, dimmx);
                const float c26 = cell_coeff_ARTM_BL (cc26, z, x, y, dimmz, dimmx);
                const float c33 = cell_coeff_BL      (cc33, z, x, y, dimmz, dimmx);
                const float c34 = cell_coeff_ARTM_BL (cc34, z, x, y, dimmz, dimmx);
                const float c35 = cell_coeff_ARTM_BL (cc35, z, x, y, dimmz, dimmx);
                const float c36 = cell_coeff_ARTM_BL (cc36, z, x, y, dimmz, dimmx);
                const float c44 = cell_coeff_BL      (cc44, z, x, y, dimmz, dimmx);
                const float c45 = cell_coeff_ARTM_BL (cc45, z, x, y, dimmz, dimmx);
                const float c46 = cell_coeff_ARTM_BL (cc46, z, x, y, dimmz, dimmx);
                const float c55 = cell_coeff_BL      (cc55, z, x, y, dimmz, dimmx);
                const float c56 = cell_coeff_ARTM_BL (cc56, z, x, y, dimmz, dimmx);
                const float c66 = cell_coeff_BL      (cc66, z, x, y, dimmz, dimmx);
                
                const float u_x = stencil_X (SX, vxu, dxi, z, x, y, dimmz, dimmx);
                const float v_x = stencil_X (SX, vxv, dxi, z, x, y, dimmz, dimmx);
                const float w_x = stencil_X (SX, vxw, dxi, z, x, y, dimmz, dimmx);
                
                const float u_y = stencil_Y (SY, vyu, dyi, z, x, y, dimmz, dimmx);
                const float v_y = stencil_Y (SY, vyv, dyi, z, x, y, dimmz, dimmx);
                const float w_y = stencil_Y (SY, vyw, dyi, z, x, y, dimmz, dimmx);
                
                const float u_z = stencil_Z (SZ, vzu, dzi, z, x, y, dimmz, dimmx);
                const float v_z = stencil_Z (SZ, vzv, dzi, z, x, y, dimmz, dimmx);
                const float w_z = stencil_Z (SZ, vzw, dzi, z, x, y, dimmz, dimmx);
                
                stress_update (sxxptr,c11,c12,c13,c14,c15,c16,z,x,y,dt,u_x,u_y,u_z,v_x,v_y,v_z,w_x,w_y,w_z,dimmz,dimmx );
                stress_update (syyptr,c12,c22,c23,c24,c25,c26,z,x,y,dt,u_x,u_y,u_z,v_x,v_y,v_z,w_x,w_y,w_z,dimmz,dimmx );
                stress_update (szzptr,c13,c23,c33,c34,c35,c36,z,x,y,dt,u_x,u_y,u_z,v_x,v_y,v_z,w_x,w_y,w_z,dimmz,dimmx );
                stress_update (syzptr,c14,c24,c34,c44,c45,c46,z,x,y,dt,u_x,u_y,u_z,v_x,v_y,v_z,w_x,w_y,w_z,dimmz,dimmx );
                stress_update (sxzptr,c15,c25,c35,c45,c55,c56,z,x,y,dt,u_x,u_y,u_z,v_x,v_y,v_z,w_x,w_y,w_z,dimmz,dimmx );
                stress_update (sxyptr,c16,c26,c36,c46,c56,c66,z,x,y,dt,u_x,u_y,u_z,v_x,v_y,v_z,w_x,w_y,w_z,dimmz,dimmx );
            }
        }
    }
}

extern "C"
void compute_component_scell_BL_cuda ( float* sxxptr,
                                       float* syyptr,
                                       float* szzptr,
                                       float* syzptr,
                                       float* sxzptr,
                                       float* sxyptr,
                                 const float* vxu,
                                 const float* vxv,
                                 const float* vxw,
                                 const float* vyu,
                                 const float* vyv,
                                 const float* vyw,
                                 const float* vzu,
                                 const float* vzv,
                                 const float* vzw,
                                 const float* cc11,
                                 const float* cc12,
                                 const float* cc13,
                                 const float* cc14,
                                 const float* cc15,
                                 const float* cc16,
                                 const float* cc22,
                                 const float* cc23,
                                 const float* cc24,
                                 const float* cc25,
                                 const float* cc26,
                                 const float* cc33,
                                 const float* cc34,
                                 const float* cc35,
                                 const float* cc36,
                                 const float* cc44,
                                 const float* cc45,
                                 const float* cc46,
                                 const float* cc55,
                                 const float* cc56,
                                 const float* cc66,
                                 const float  dt,
                                 const float  dzi,
                                 const float  dxi,
                                 const float  dyi,
                                 const int    nz0,
                                 const int    nzf,
                                 const int    nx0,
                                 const int    nxf,
                                 const int    ny0,
                                 const int    nyf,
                                 const int    SZ,
                                 const int    SX,
                                 const int    SY,
                                 const int    dimmz,
                                 const int    dimmx,
                                 void*        stream)
{
    const int block_dim_x = 32;
    const int block_dim_y = 4;


    dim3 grid_dim( ((nzf-nz0) + block_dim_x-1)/block_dim_x,
                   ((nxf-nx0) + block_dim_y-1)/block_dim_y,
                      1 );
    dim3 block_dim(block_dim_x, block_dim_y, 1);

    cudaStream_t s = (cudaStream_t) stream;

    compute_component_scell_BL_cuda_k<<<grid_dim, block_dim, 0, s>>>
        (sxxptr, syyptr, szzptr, syzptr, sxzptr, sxyptr,
         vxu, vxv, vxw, vyu, vyv, vyw, vzu, vzv, vzw,
         cc11, cc12, cc13, cc14, cc15, cc16,
         cc22, cc23, cc24, cc25, cc26,
         cc33, cc34, cc35, cc36,
         cc44, cc45, cc45,
         cc55, cc56,
         cc66,
         dt, dzi, dxi, dyi, nz0, nzf, nx0, nxf, ny0, nyf, SZ, SX, SY, dimmz, dimmx);
};



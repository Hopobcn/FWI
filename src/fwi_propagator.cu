#include <cuda_runtime.h>
#include <cstdio>

#include "fwi_propagator.cuh"

#define C0 1.2f
#define C1 1.4f
#define C2 1.6f
#define C3 1.8f

cudaTextureObject_t create_texture_obj(const float* d_data,
                                       cudaArray_t& a_data,
                                       const dim_t dim,
                                       cudaTextureFilterMode filterMode = cudaFilterModeLinear,
                                       cudaTextureAddressMode addressMode = cudaAddressModeClamp);

__host__
cudaTextureObject_t create_texture_obj(const float* d_data,
                                       cudaArray_t& a_data,
                                       const dim_t dim,
                                       cudaTextureFilterMode filterMode,
                                       cudaTextureAddressMode addressMode)
{
    // create channel descriptor
    cudaChannelFormatDesc channelDesc = cudaCreateChannelDesc<float>();

    // create cuda array
    cudaExtent extent_elems = make_cudaExtent( dim.zsize              , dim.xsize, dim.ysize );
    CUDA_CHECK( cudaMalloc3DArray(&a_data, &channelDesc, extent_elems) );
    
    // 2.b copy data to container
    cudaMemcpy3DParms params = {0};
    params.srcPtr   = make_cudaPitchedPtr((void*)d_data, dim.pitch*sizeof(float), dim.xsize, dim.ysize);
    params.dstArray = a_data;
    params.extent   = extent_elems;   // dimensions of the transfered area in elements
    params.kind     = cudaMemcpyDeviceToDevice;

    CUDA_CHECK( cudaMemcpy3D(&params) );

    // specify texture contents
    struct cudaResourceDesc resDesc;
    memset(&resDesc, 0, sizeof(resDesc));
    resDesc.resType = cudaResourceTypeArray;
    resDesc.res.array.array = a_data;

    // texture obj parameters
    struct cudaTextureDesc texDesc;
    memset(&texDesc, 0, sizeof(texDesc));
    texDesc.addressMode[0] = addressMode;
    texDesc.addressMode[1] = addressMode;
    texDesc.addressMode[2] = addressMode;
    texDesc.filterMode     = filterMode;
    texDesc.readMode       = cudaReadModeElementType;
    texDesc.normalizedCoords = false;

    // create texture object
    cudaTextureObject_t texObj = 0;
    CUDA_CHECK( cudaCreateTextureObject(&texObj, &resDesc, &texDesc, NULL) );

    return texObj;
}

void destroy_texture_obj(cudaTextureObject_t obj, cudaArray_t a_data)
{
    CUDA_CHECK( cudaDestroyTextureObject(obj) );
    CUDA_CHECK( cudaFreeArray(a_data) );
}

__device__ inline
int IDX (const int   z, 
         const int   x, 
         const int   y, 
         const dim_t dim)
{
    return (y*dim.xsize + x) * dim.pitch + z;
};

__device__ inline
float stencil_Z (const int   off,
                 const float* __restrict__ ptr,
                 const float dzi,
                 const int   z,
                 const int   x,
                 const int   y,
                 const dim_t dim)
{
    return  ((C0 * ( ptr[IDX(z  +off,x,y,dim)] - ptr[IDX(z-1+off,x,y,dim)]) +
              C1 * ( ptr[IDX(z+1+off,x,y,dim)] - ptr[IDX(z-2+off,x,y,dim)]) +
              C2 * ( ptr[IDX(z+2+off,x,y,dim)] - ptr[IDX(z-3+off,x,y,dim)]) +
              C3 * ( ptr[IDX(z+3+off,x,y,dim)] - ptr[IDX(z-4+off,x,y,dim)])) * dzi );
};

template <const int HALO,
          const int BDIMX>
__device__ inline
float stencil_Z_shfl (const int   off,
                      const float* __restrict__ ptr_gmem,
                      const float dzi,
                      const int   z,
                      const int   x,
                      const int   y,
                      const dim_t dim)
{
    float current = (z+off < dim.zsize && x < dim.xsize) ? ptr_gmem[IDX(z+off,x,y,dim)] : 0.0f;
    float right3  = __shfl_down(current, 3);
    float right2  = __shfl_down(current, 2);
    float right1  = __shfl_down(current, 1);
    float left1   = __shfl_up(current, 1);
    float left2   = __shfl_up(current, 2);
    float left3   = __shfl_up(current, 3);
    float left4   = __shfl_up(current, 4);

    /* For threads without neighbors: */
    if (threadIdx.x < 1 /* 1 */) left1 = (z+off-1 < dim.zsize && x < dim.xsize) ? ptr_gmem[IDX(z+off-1,x,y,dim)] : 0.0f;
    if (threadIdx.x < 2 /* 2 */) left2 = (z+off-2 < dim.zsize && x < dim.xsize) ? ptr_gmem[IDX(z+off-2,x,y,dim)] : 0.0f;
    if (threadIdx.x < 3 /* 3 */) left3 = (z+off-3 < dim.zsize && x < dim.xsize) ? ptr_gmem[IDX(z+off-3,x,y,dim)] : 0.0f;
    if (threadIdx.x < 4 /* 4 */) left4 = (z+off-4 < dim.zsize && x < dim.xsize) ? ptr_gmem[IDX(z+off-4,x,y,dim)] : 0.0f;

    if (threadIdx.x >= BDIMX-1 /* 1 */) right1 = (z+off+1 < dim.zsize && x < dim.xsize) ? ptr_gmem[IDX(z+off+1,x,y,dim)] : 0.0f;
    if (threadIdx.x >= BDIMX-2 /* 2 */) right2 = (z+off+1 < dim.zsize && x < dim.xsize) ? ptr_gmem[IDX(z+off+2,x,y,dim)] : 0.0f;
    if (threadIdx.x >= BDIMX-3 /* 3 */) right3 = (z+off+1 < dim.zsize && x < dim.xsize) ? ptr_gmem[IDX(z+off+3,x,y,dim)] : 0.0f;

    return  ((C0 * ( current - left1) +
              C1 * ( right1  - left2) +
              C2 * ( right2  - left3) +
              C3 * ( right3  - left4)) * dzi );
};

__device__ inline
float stencil_X(const int   off,
                const float* __restrict__ ptr,
                const float dxi,
                const int   z,
                const int   x,
                const int   y,
                const dim_t dim)
{
    return ((C0 * ( ptr[IDX(z,x  +off,y,dim)] - ptr[IDX(z,x-1+off,y,dim)]) +
             C1 * ( ptr[IDX(z,x+1+off,y,dim)] - ptr[IDX(z,x-2+off,y,dim)]) +
             C2 * ( ptr[IDX(z,x+2+off,y,dim)] - ptr[IDX(z,x-3+off,y,dim)]) +
             C3 * ( ptr[IDX(z,x+3+off,y,dim)] - ptr[IDX(z,x-4+off,y,dim)])) * dxi );
};



template <const int NX,
          const int NY,
          const int HALO,
          const int BDIMY>
__device__ inline
float stencil_X_smem(const int   off,
                           float ptr_smem[NY][NX],
                     const float* __restrict__ ptr_gmem,
                     const float dxi,
                     const int   z,
                     const int   x,
                     const int   y,
                     const dim_t dim)
{
    //__syncthreads();
    const int tx = threadIdx.x;
    const int ty = threadIdx.y+HALO;
    ///////////// intra-block communication///////////////////
    ptr_smem[ty][tx] = ((z < dim.zsize && x+off < dim.xsize) ? ptr_gmem[IDX(z,x+off,y,dim)] : 0.0f);
    if (threadIdx.y < HALO)
    {
        ptr_smem[ty-HALO ][tx] = ((z < dim.zsize && x+off-HALO  < dim.xsize) ? ptr_gmem[IDX(z,x+off-HALO, y,dim)] : 0.0f);
        ptr_smem[ty+BDIMY][tx] = ((z < dim.zsize && x+off+BDIMY < dim.xsize) ? ptr_gmem[IDX(z,x+off+BDIMY,y,dim)] : 0.0f);
    }
    __syncthreads();

    return ((C0 * ( ptr_smem[ty  ][tx] - ptr_smem[ty-1][tx] ) +
             C1 * ( ptr_smem[ty+1][tx] - ptr_smem[ty-2][tx] ) +
             C2 * ( ptr_smem[ty+2][tx] - ptr_smem[ty-3][tx] ) +
             C3 * ( ptr_smem[ty+3][tx] - ptr_smem[ty-4][tx] )) * dxi );
};

__device__ inline
float stencil_Y(const int   off,
                const float* __restrict__ ptr,
                const float dyi,
                const int   z,
                const int   x,
                const int   y,
                const dim_t dim)
{
    return ((C0 * ( ptr[IDX(z,x,y  +off,dim)] - ptr[IDX(z,x,y-1+off,dim)]) +
             C1 * ( ptr[IDX(z,x,y+1+off,dim)] - ptr[IDX(z,x,y-2+off,dim)]) +
             C2 * ( ptr[IDX(z,x,y+2+off,dim)] - ptr[IDX(z,x,y-3+off,dim)]) +
             C3 * ( ptr[IDX(z,x,y+3+off,dim)] - ptr[IDX(z,x,y-4+off,dim)])) * dyi );
};

/* -------------------------------------------------------------------- */
/*                     KERNELS FOR VELOCITY                             */
/* -------------------------------------------------------------------- */

__device__ inline
float rho_BL (const float* __restrict__ rho,
              const int   z,
              const int   x,
              const int   y,
              const dim_t dim)
{
    return (2.0f / (rho[IDX(z,x,y,dim)] + rho[IDX(z+1,x,y,dim)]));
};

template <const int BDIMX>
__device__ inline
float rho_BL_shfl (const float* __restrict__ ptr_gmem,
                   const int   z,
                   const int   x,
                   const int   y,
                   const dim_t dim)
{
    float current = (z<dim.zsize && x<dim.xsize) ? ptr_gmem[IDX(z,x,y,dim)] : 0.0f;
    float right1  = __shfl_down(current, 1);

    /* For threads without neighbors: */
    if (threadIdx.x >= BDIMX-1 /* 1 */) right1 = (z+1<dim.zsize && x<dim.xsize) ? ptr_gmem[IDX(z+1,x,y,dim)] : 0.0f;

    return (2.0f / (current + right1));
};

__device__ inline
float rho_TR (const float* __restrict__ rho,
              const int   z,
              const int   x,
              const int   y,
              const dim_t dim)
{
    return (2.0f / (rho[IDX(z,x,y,dim)] + rho[IDX(z,x+1,y,dim)]));
};

template <const int BDIMX,
          const int BDIMY>
__device__ inline
float rho_TR_smem (      float rho_smem[BDIMY+1][BDIMX],
                   const float* __restrict__ rho_gmem,
                   const int   z,
                   const int   x,
                   const int   y,
                   const dim_t dim)
{
    const int tx = threadIdx.x;
    const int ty = threadIdx.y;

    rho_smem[ty][tx] = (z<dim.zsize && x<dim.xsize) ? rho_gmem[IDX(z,x,y,dim)] : 0.0f;

    /* For threads without neighbors: */
    if (ty < 1)
        rho_smem[ty+BDIMY][tx] = (z<dim.zsize && x+BDIMY<dim.xsize) ? rho_gmem[IDX(z,x+BDIMY,y,dim)] : 0.0f;

    __syncthreads();
    
    return (2.0f / (rho_smem[ty][tx] + rho_smem[ty+1][tx]));
};


__device__ inline
float rho_BR (const float* __restrict__ rho,
              const int   z,
              const int   x,
              const int   y,
              const dim_t dim)
{
    return ( 8.0f/ ( rho[IDX(z  ,x  ,y  ,dim)] +
                     rho[IDX(z+1,x  ,y  ,dim)] +
                     rho[IDX(z  ,x+1,y  ,dim)] +
                     rho[IDX(z  ,x  ,y+1,dim)] +
                     rho[IDX(z  ,x+1,y+1,dim)] +
                     rho[IDX(z+1,x+1,y  ,dim)] +
                     rho[IDX(z+1,x  ,y+1,dim)] +
                     rho[IDX(z+1,x+1,y+1,dim)]) );
};

__device__ inline
float rho_BR (cudaTextureObject_t tex,
              const float z,
              const float x,
              const float y)
{
    return 1.0f/tex3D<float>(tex, z+0.5f, x+0.5f, y+0.5f);
};

template <const int BDIMX,
          const int BDIMY>
__device__ inline
float rho_BR_smem (      float rho_smem[BDIMY+1][BDIMX+1],
                   const float* __restrict__ rho_gmem,
                   const float rho_current,
                   const float rho_front1,
                   const int   z,
                   const int   x,
                   const int   y,
                   const dim_t dim)
{
    const int tx = threadIdx.x;
    const int ty = threadIdx.y;

    rho_smem[ty][tx] = (z<dim.zsize && x<dim.xsize) ? rho_gmem[IDX(z,x,y,dim)] : 0.0f;
    
    /* For threads without neighbors: */
    if (ty < 1)
        rho_smem[ty+BDIMY][tx      ] = (z<dim.zsize && x+BDIMY<dim.xsize) ? rho_gmem[IDX(z,x+BDIMY,y,dim)] : 0.0f;
    if (tx < 1)
        rho_smem[ty      ][tx+BDIMX] = (z+BDIMX<dim.zsize && x<dim.xsize) ? rho_gmem[IDX(z+BDIMX,x,y,dim)] : 0.0f;
    
    __syncthreads();
   
    return (8.0f/ ( rho_current                    +
                    rho_smem[ty  ][tx+1]           +
                    rho_smem[ty+1][tx  ]           +
                    rho_front1                     +
                    rho_gmem[IDX(z  ,x+1,y+1,dim)] +
                    rho_gmem[IDX(z+1,x+1,y  ,dim)] +
                    rho_gmem[IDX(z+1,x  ,y+1,dim)] +
                    rho_gmem[IDX(z+1,x+1,y+1,dim)]) );
};

__device__ inline
float rho_TL (const float* __restrict__ rho,
              const int   z,
              const int   x,
              const int   y,
              const dim_t dim)
{
    return (2.0f / (rho[IDX(z,x,y,dim)] + rho[IDX(z,x,y+1,dim)]));
};

#ifdef VCELL_TL
template <const int HALO  = 4, 
          const int BDIMX = 32,
          const int BDIMY = 4,
          const int NX    = BDIMX+2*HALO,
          const int NY    = BDIMY+2*HALO>
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
                                   const dim_t           dim)
{
    int z = blockIdx.x * blockDim.x + threadIdx.x + nz0; 
    int x = blockIdx.y * blockDim.y + threadIdx.y + nx0; 
    
    // WARNING: We can't predicate threads that fall outside of the [nz0:nzf][nx0:nxf] range because
    //          we use COLLECTIVE operations like SHUFFLE & SHARED.
    //          PREVENT incorrect GMEM memory access by checking boundaries at every access

    __shared__ float sx_smem[NY][NX];
    float sy_front1, sy_front2, sy_front3;
    float sy_back1, sy_back2, sy_back3, sy_back4;
    float sy_current;

    sy_back3   = (z<dim.zsize && x<dim.xsize) ? syptr[IDX(z,x,ny0-HALO+0+SY,dim)] : 0.0f;
    sy_back2   = (z<dim.zsize && x<dim.xsize) ? syptr[IDX(z,x,ny0-HALO+1+SY,dim)] : 0.0f;
    sy_back1   = (z<dim.zsize && x<dim.xsize) ? syptr[IDX(z,x,ny0-HALO+2+SY,dim)] : 0.0f;
    sy_current = (z<dim.zsize && x<dim.xsize) ? syptr[IDX(z,x,ny0-HALO+3+SY,dim)] : 0.0f;
    sy_front1  = (z<dim.zsize && x<dim.xsize) ? syptr[IDX(z,x,ny0-HALO+4+SY,dim)] : 0.0f;
    sy_front2  = (z<dim.zsize && x<dim.xsize) ? syptr[IDX(z,x,ny0-HALO+5+SY,dim)] : 0.0f;
    sy_front3  = (z<dim.zsize && x<dim.xsize) ? syptr[IDX(z,x,ny0-HALO+6+SY,dim)] : 0.0f;

    float rho_current, rho_front1;
    rho_front1 = (z<dim.zsize && x<dim.xsize) ? rho[IDX(z,x,ny0,dim)] : 0.0f;

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
        sy_front3  = (z<dim.zsize && x<dim.xsize) ? syptr[IDX(z,x,y+SY+HALO-1,dim)] : 0.0f;
        ///////////////////////
        rho_current = rho_front1;
        rho_front1  = (z<dim.zsize && x<dim.xsize) ? rho[IDX(z,x,y+1,dim)] : 0.0f;
        //////////////////////////////////////////////////////////

        const float lrho = (2.0f / (rho_current + rho_front1));

        const float stz = stencil_Z_shfl <HALO,BDIMX> (SZ, szptr, dzi, z, x, y, dim);
        
        const float stx = stencil_X_smem <NX,NY,HALO,BDIMY> (SX, sx_smem, sxptr, dxi, z, x, y, dim);

        const float sty = ((C0 * ( sy_current - sy_back1 ) +
                            C1 * ( sy_front1  - sy_back2 ) +
                            C2 * ( sy_front2  - sy_back3 ) +
                            C3 * ( sy_front3  - sy_back4 )) * dyi );
    
        if (z < nzf && x < nxf)
        {
            vptr[IDX(z,x,y,dim)] += (stx  + sty  + stz) * dt * lrho;
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
                                   const dim_t           dim)
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
                const float lrho = rho_TL(rho, z, x, y, dim);
                
                const float stx  = stencil_X( SX, sxptr, dxi, z, x, y, dim);
                const float sty  = stencil_Y( SY, syptr, dyi, z, x, y, dim);
                const float stz  = stencil_Z( SZ, szptr, dzi, z, x, y, dim);
                
                vptr[IDX(z,x,y,dim)] += (stx  + sty  + stz) * dt * lrho;
            }
        }
    }
}
#endif

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
                                 const dim_t  dim,
                                 void*        stream)
{
    const int block_dim_x = 32;
    const int block_dim_y = 4;


    dim3 grid_dim( ((nzf-nz0) + block_dim_x-1)/block_dim_x,
                   ((nxf-nx0) + block_dim_y-1)/block_dim_y,
                      1 );
    dim3 block_dim(block_dim_x, block_dim_y, 1);

    cudaStream_t s = (cudaStream_t) stream;

#ifdef VCELL_TL
    compute_component_vcell_TL_cuda_k<4,block_dim_x,block_dim_y><<<grid_dim, block_dim, 0, s>>>
        (vptr, szptr, sxptr, syptr, rho, dt, dzi, dxi, dyi, 
         nz0, nzf, nx0, nxf, ny0, nyf, SZ, SX, SY, dim);
#else
    compute_component_vcell_TL_cuda_k<<<grid_dim, block_dim, 0, s>>>
        (vptr, szptr, sxptr, syptr, rho, dt, dzi, dxi, dyi, 
         nz0, nzf, nx0, nxf, ny0, nyf, SZ, SX, SY, dim);
#endif
    CUDA_CHECK(cudaGetLastError());
};


#ifdef VCELL_TR
template <const int HALO  = 4, 
          const int BDIMX = 32,
          const int BDIMY = 4,
          const int NX    = BDIMX+2*HALO,
          const int NY    = BDIMY+2*HALO>
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
                                   const dim_t           dim)
{
    int z = blockIdx.x * blockDim.x + threadIdx.x + nz0; 
    int x = blockIdx.y * blockDim.y + threadIdx.y + nx0; 
    
    // WARNING: We can't predicate threads that fall outside of the [nz0:nzf][nx0:nxf] range because
    //          we use COLLECTIVE operations like SHUFFLE & SHARED.
    //          PREVENT incorrect GMEM memory access by checking boundaries at every access

    __shared__ float sx_smem[NY][NX];
    float sy_front1, sy_front2, sy_front3;
    float sy_back1, sy_back2, sy_back3, sy_back4;
    float sy_current;

    sy_back3   = (z<dim.zsize && x<dim.xsize) ? syptr[IDX(z,x,ny0-HALO+0+SY,dim)] : 0.0f;
    sy_back2   = (z<dim.zsize && x<dim.xsize) ? syptr[IDX(z,x,ny0-HALO+1+SY,dim)] : 0.0f;
    sy_back1   = (z<dim.zsize && x<dim.xsize) ? syptr[IDX(z,x,ny0-HALO+2+SY,dim)] : 0.0f;
    sy_current = (z<dim.zsize && x<dim.xsize) ? syptr[IDX(z,x,ny0-HALO+3+SY,dim)] : 0.0f;
    sy_front1  = (z<dim.zsize && x<dim.xsize) ? syptr[IDX(z,x,ny0-HALO+4+SY,dim)] : 0.0f;
    sy_front2  = (z<dim.zsize && x<dim.xsize) ? syptr[IDX(z,x,ny0-HALO+5+SY,dim)] : 0.0f;
    sy_front3  = (z<dim.zsize && x<dim.xsize) ? syptr[IDX(z,x,ny0-HALO+6+SY,dim)] : 0.0f;

    __shared__ float rho_smem[BDIMY+1][BDIMX];

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
        sy_front3  = (z<dim.zsize && x<dim.xsize) ? syptr[IDX(z,x,y+SY+HALO-1,dim)] : 0.0f;
        //////////////////////////////////////////////////////////

        const float lrho = rho_TR_smem <BDIMX,BDIMY> (rho_smem, rho, z, x, y, dim);

        const float stz = stencil_Z_shfl <HALO,BDIMX> (SZ, szptr, dzi, z, x, y, dim);
        
        const float stx = stencil_X_smem <NX,NY,HALO,BDIMY> (SX, sx_smem, sxptr, dxi, z, x, y, dim);

        const float sty = ((C0 * ( sy_current - sy_back1 ) +
                            C1 * ( sy_front1  - sy_back2 ) +
                            C2 * ( sy_front2  - sy_back3 ) +
                            C3 * ( sy_front3  - sy_back4 )) * dyi );
        
        if (z < nzf && x < nxf)
        {
            vptr[IDX(z,x,y,dim)] += (stx  + sty  + stz) * dt * lrho;
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
                                   const dim_t           dim)
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
                const float lrho = rho_TR(rho, z, x, y, dim);
                
                const float stx  = stencil_X( SX, sxptr, dxi, z, x, y, dim);
                const float sty  = stencil_Y( SY, syptr, dyi, z, x, y, dim);
                const float stz  = stencil_Z( SZ, szptr, dzi, z, x, y, dim);
                
                vptr[IDX(z,x,y,dim)] += (stx  + sty  + stz) * dt * lrho;
            }
        }
    }
}
#endif

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
                                 const dim_t  dim,
                                 void*        stream)
{
    const int block_dim_x = 32;
    const int block_dim_y = 4;


    dim3 grid_dim( ((nzf-nz0) + block_dim_x-1)/block_dim_x,
                   ((nxf-nx0) + block_dim_y-1)/block_dim_y,
                      1 );
    dim3 block_dim(block_dim_x, block_dim_y, 1);

    cudaStream_t s = (cudaStream_t) stream;

#ifdef VCELL_TR
    compute_component_vcell_TR_cuda_k<4,block_dim_x,block_dim_y><<<grid_dim, block_dim, 0, s>>>
        (vptr, szptr, sxptr, syptr, rho, dt, dzi, dxi, dyi, 
         nz0, nzf, nx0, nxf, ny0, nyf, SZ, SX, SY, dim);
#else
    compute_component_vcell_TR_cuda_k<<<grid_dim, block_dim, 0, s>>>
        (vptr, szptr, sxptr, syptr, rho, dt, dzi, dxi, dyi, 
         nz0, nzf, nx0, nxf, ny0, nyf, SZ, SX, SY, dim);
#endif 
    CUDA_CHECK(cudaGetLastError());
};


#ifdef VCELL_BR
template <const int HALO  = 4, 
          const int BDIMX = 32,
          const int BDIMY = 4,
          const int NX    = BDIMX+2*HALO,
          const int NY    = BDIMY+2*HALO>
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
                                   const dim_t           dim)
{
    int z = blockIdx.x * blockDim.x + threadIdx.x + nz0; 
    int x = blockIdx.y * blockDim.y + threadIdx.y + nx0; 
    
    __shared__ float sx_smem[NY][NX];
    float sy_front1, sy_front2, sy_front3;
    float sy_back1, sy_back2, sy_back3, sy_back4;
    float sy_current;

    sy_back3   = (z<dim.zsize && x<dim.xsize) ? syptr[IDX(z,x,ny0-HALO+0+SY,dim)] : 0.0f;
    sy_back2   = (z<dim.zsize && x<dim.xsize) ? syptr[IDX(z,x,ny0-HALO+1+SY,dim)] : 0.0f;
    sy_back1   = (z<dim.zsize && x<dim.xsize) ? syptr[IDX(z,x,ny0-HALO+2+SY,dim)] : 0.0f;
    sy_current = (z<dim.zsize && x<dim.xsize) ? syptr[IDX(z,x,ny0-HALO+3+SY,dim)] : 0.0f;
    sy_front1  = (z<dim.zsize && x<dim.xsize) ? syptr[IDX(z,x,ny0-HALO+4+SY,dim)] : 0.0f;
    sy_front2  = (z<dim.zsize && x<dim.xsize) ? syptr[IDX(z,x,ny0-HALO+5+SY,dim)] : 0.0f;
    sy_front3  = (z<dim.zsize && x<dim.xsize) ? syptr[IDX(z,x,ny0-HALO+6+SY,dim)] : 0.0f;

    __shared__ float rho_smem[BDIMY+1][BDIMX+1];
    float rho_current, rho_front1;
    rho_front1 = (z<dim.zsize && x<dim.xsize) ? rho[IDX(z,x,ny0,dim)] : 0.0f;

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
        sy_front3  = (z<dim.zsize && x<dim.xsize) ? syptr[IDX(z,x,y+SY+HALO-1,dim)] : 0.0f;
        ///////////////////////
        rho_current = rho_front1;
        rho_front1  = (z<dim.zsize && x<dim.xsize) ? rho[IDX(z,x,y+1,dim)] : 0.0f;
        //////////////////////////////////////////////////////////

        const float lrho = rho_BR_smem <BDIMX,BDIMY> (rho_smem, rho, rho_current, rho_front1, z, x, y, dim);
        
        const float stz = stencil_Z_shfl <HALO,BDIMX> (SZ, szptr, dzi, z, x, y, dim);

        const float stx = stencil_X_smem <NX,NY,HALO,BDIMY> (SX, sx_smem, sxptr, dxi, z, x, y, dim);

        const float sty = ((C0 * ( sy_current - sy_back1 ) +
                            C1 * ( sy_front1  - sy_back2 ) +
                            C2 * ( sy_front2  - sy_back3 ) +
                            C3 * ( sy_front3  - sy_back4 )) * dyi );

        if (z < nzf && x < nxf)
        {
            vptr[IDX(z,x,y,dim)] += (stx  + sty  + stz) * dt * lrho;
        }
    }
}
#elif defined(VCELL_BR_TEXTURE)
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
                                   const dim_t           dim,
                                   cudaTextureObject_t   tex)
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
                const float lrho    = rho_BR(tex, z, x, y);

                const float stx  = stencil_X( SX, sxptr, dxi, z, x, y, dim);
                const float sty  = stencil_Y( SY, syptr, dyi, z, x, y, dim);
                const float stz  = stencil_Z( SZ, szptr, dzi, z, x, y, dim);
                
                vptr[IDX(z,x,y,dim)] += (stx  + sty  + stz) * dt * lrho;
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
                                   const dim_t           dim)
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
                const float lrho = rho_BR(rho, z, x, y, dim);
                
                const float stx  = stencil_X( SX, sxptr, dxi, z, x, y, dim);
                const float sty  = stencil_Y( SY, syptr, dyi, z, x, y, dim);
                const float stz  = stencil_Z( SZ, szptr, dzi, z, x, y, dim);
                
                vptr[IDX(z,x,y,dim)] += (stx  + sty  + stz) * dt * lrho;
            }
        }
    }
}
#endif

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
                                 const dim_t  dim,
                                 void*        stream)
{
    const int block_dim_x = 32;
    const int block_dim_y = 4;


    dim3 grid_dim( ((nzf-nz0) + block_dim_x-1)/block_dim_x,
                   ((nxf-nx0) + block_dim_y-1)/block_dim_y,
                      1 );
    dim3 block_dim(block_dim_x, block_dim_y, 1);

    cudaStream_t s = (cudaStream_t) stream;

#ifdef VCELL_BR_TEXTURE
    cudaArray_t a_data;
    cudaTextureObject_t tex = create_texture_obj(rho, a_data, dim);
#endif

#ifdef VCELL_BR
    compute_component_vcell_BR_cuda_k<4,block_dim_x,block_dim_y><<<grid_dim, block_dim, 0, s>>>
        (vptr, szptr, sxptr, syptr, rho, dt, dzi, dxi, dyi, 
         nz0, nzf, nx0, nxf, ny0, nyf, SZ, SX, SY, dim);
#else
#ifdef VCELL_BR_TEXTURE
    compute_component_vcell_BR_cuda_k<<<grid_dim, block_dim, 0, s>>>
        (vptr, szptr, sxptr, syptr, rho, dt, dzi, dxi, dyi, 
         nz0, nzf, nx0, nxf, ny0, nyf, SZ, SX, SY, dim, tex);
#else
    compute_component_vcell_BR_cuda_k<<<grid_dim, block_dim, 0, s>>>
        (vptr, szptr, sxptr, syptr, rho, dt, dzi, dxi, dyi, 
         nz0, nzf, nx0, nxf, ny0, nyf, SZ, SX, SY, dim);
#endif
#endif 
    CUDA_CHECK(cudaGetLastError());

#ifdef VCELL_BR_TEXTURE
    destroy_texture_obj(tex, a_data);
#endif
};

#ifdef VCELL_BL
template <const int HALO  = 4, 
          const int BDIMX = 32,
          const int BDIMY = 4,
          const int NX    = BDIMX+2*HALO,
          const int NY    = BDIMY+2*HALO>
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
                                   const dim_t           dim)
{
    int z = blockIdx.x * blockDim.x + threadIdx.x + nz0; 
    int x = blockIdx.y * blockDim.y + threadIdx.y + nx0; 
    
    __shared__ float sx_smem[NY][NX];
    float sy_front1, sy_front2, sy_front3;
    float sy_back1, sy_back2, sy_back3, sy_back4;
    float sy_current;

    sy_back3   = (z<dim.zsize && x<dim.xsize) ? syptr[IDX(z,x,ny0-HALO+0+SY,dim)] : 0.0f;
    sy_back2   = (z<dim.zsize && x<dim.xsize) ? syptr[IDX(z,x,ny0-HALO+1+SY,dim)] : 0.0f;
    sy_back1   = (z<dim.zsize && x<dim.xsize) ? syptr[IDX(z,x,ny0-HALO+2+SY,dim)] : 0.0f;
    sy_current = (z<dim.zsize && x<dim.xsize) ? syptr[IDX(z,x,ny0-HALO+3+SY,dim)] : 0.0f;
    sy_front1  = (z<dim.zsize && x<dim.xsize) ? syptr[IDX(z,x,ny0-HALO+4+SY,dim)] : 0.0f;
    sy_front2  = (z<dim.zsize && x<dim.xsize) ? syptr[IDX(z,x,ny0-HALO+5+SY,dim)] : 0.0f;
    sy_front3  = (z<dim.zsize && x<dim.xsize) ? syptr[IDX(z,x,ny0-HALO+6+SY,dim)] : 0.0f;

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
        sy_front3  = (z<dim.zsize && x<dim.xsize) ? syptr[IDX(z,x,y+SY+HALO-1,dim)] : 0.0f;
        //////////////////////////////////////////////////////////
        
        const float lrho = rho_BL_shfl <BDIMX> (rho, z, x, y, dim);

        const float stz = stencil_Z_shfl <HALO,BDIMX> (SZ, szptr, dzi, z, x, y, dim);

        const float stx = stencil_X_smem <NX,NY,HALO,BDIMY> (SX, sx_smem, sxptr, dxi, z, x, y, dim);

        const float sty = ((C0 * ( sy_current - sy_back1 ) +
                            C1 * ( sy_front1  - sy_back2 ) +
                            C2 * ( sy_front2  - sy_back3 ) +
                            C3 * ( sy_front3  - sy_back4 )) * dyi );

        if (z < nzf && x < nxf)
        {
            vptr[IDX(z,x,y,dim)] += (stx  + sty  + stz) * dt * lrho;
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
                                   const dim_t           dim)
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
                const float lrho = rho_BL(rho, z, x, y, dim);
                
                const float stx  = stencil_X( SX, sxptr, dxi, z, x, y, dim);
                const float sty  = stencil_Y( SY, syptr, dyi, z, x, y, dim);
                const float stz  = stencil_Z( SZ, szptr, dzi, z, x, y, dim);
                
                vptr[IDX(z,x,y,dim)] += (stx  + sty  + stz) * dt * lrho;
            }
        }
    }
}
#endif

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
                                 const dim_t  dim,
                                 void*        stream)
{
    const int block_dim_x = 32;
    const int block_dim_y = 4;


    dim3 grid_dim( ((nzf-nz0) + block_dim_x-1)/block_dim_x,
                   ((nxf-nx0) + block_dim_y-1)/block_dim_y,
                      1 );
    dim3 block_dim(block_dim_x, block_dim_y, 1);

    cudaStream_t s = (cudaStream_t) stream;

#ifdef VCELL_BL
    compute_component_vcell_BL_cuda_k<4,block_dim_x,block_dim_y><<<grid_dim, block_dim, 0, s>>>
        (vptr, szptr, sxptr, syptr, rho, dt, dzi, dxi, dyi, 
         nz0, nzf, nx0, nxf, ny0, nyf, SZ, SX, SY, dim);
#else
    compute_component_vcell_BL_cuda_k<<<grid_dim, block_dim, 0, s>>>
        (vptr, szptr, sxptr, syptr, rho, dt, dzi, dxi, dyi, 
         nz0, nzf, nx0, nxf, ny0, nyf, SZ, SX, SY, dim);
#endif 
    CUDA_CHECK(cudaGetLastError());
};





/* ------------------------------------------------------------------------------ */
/*                                                                                */
/*                               CALCULO DE TENSIONES                             */
/*                                                                                */
/* ------------------------------------------------------------------------------ */

__device__ inline
void stress_update(float* __restrict__ sptr,
                   const float c1,
                   const float c2,
                   const float c3,
                   const float c4,
                   const float c5,
                   const float c6,
                   const int   z,
                   const int   x,
                   const int   y,
                   const float dt,
                   const float u_x,
                   const float u_y,
                   const float u_z,
                   const float v_x,
                   const float v_y,
                   const float v_z,
                   const float w_x,
                   const float w_y,
                   const float w_z,
                   const dim_t dim)
{
    float accum  = dt * c1 * u_x;
         accum += dt * c2 * v_y;
         accum += dt * c3 * w_z;
         accum += dt * c4 * (w_y + v_z);
         accum += dt * c5 * (w_x + u_z);
         accum += dt * c6 * (v_x + u_y);
    sptr[IDX(z,x,y,dim)] += accum;
};

__device__ inline
float cell_coeff_BR (const float* __restrict__ ptr, 
                     const int   z,
                     const int   x,
                     const int   y,
                     const dim_t dim)
{
    return ( 1.0f / ( 2.5f  *(ptr[IDX(z  , x  ,y,dim)] +
                              ptr[IDX(z  , x+1,y,dim)] +
                              ptr[IDX(z+1, x  ,y,dim)] +
                              ptr[IDX(z+1, x+1,y,dim)])) );
};

template <const int NX,
          const int NY,
          const int BDIMX,
          const int BDIMY>
__device__ inline
float cell_coeff_BR_smem (      float ptr_smem[NY][NX],
                          const float* __restrict__ ptr_gmem, 
                          const int   z,
                          const int   x,
                          const int   y,
                          const dim_t dim)
{
    const int tx = threadIdx.x;
    const int ty = threadIdx.y;
    ///////////// intra-block communication///////////////////
    ptr_smem[ty][tx] = (z < dim.zsize && x < dim.xsize) ? ptr_gmem[IDX(z,x,y,dim)] : 0.0f;
    if (tx < 1) ptr_smem[ty][tx+BDIMX] = (z+BDIMX < dim.zsize && x < dim.xsize) ? ptr_gmem[IDX(z+BDIMX,x,y,dim)] : 0.0f;
    if (ty < 1) ptr_smem[ty+BDIMY][tx] = (z < dim.zsize && x+BDIMY < dim.xsize) ? ptr_gmem[IDX(z,x+BDIMY,y,dim)] : 0.0f;
    if (tx == 0 && ty == 0) 
                ptr_smem[ty+BDIMY][tx+BDIMX] = (z+BDIMX<dim.zsize && x+BDIMY<dim.xsize) ? ptr_gmem[IDX(z+BDIMX,x+BDIMY,y,dim)] : 0.0f;

    __syncthreads();

    return ( 1.0f / ( 2.5f  *(ptr_smem[ty  ][tx  ] +
                              ptr_smem[ty+1][tx  ] +
                              ptr_smem[ty  ][tx+1] +
                              ptr_smem[ty+1][tx+1])) );
};

__device__ inline
float cell_coeff_TL (const float* __restrict__ ptr, 
                     const int   z, 
                     const int   x, 
                     const int   y, 
                     const dim_t dim)
{
    return ( 1.0f / (ptr[IDX(z,x,y,dim)]));
};

__device__ inline
float cell_coeff_BL (const float* __restrict__ ptr, 
                     const int   z, 
                     const int   x, 
                     const int   y, 
                     const dim_t dim)
{
    return ( 1.0f / ( 2.5f *(ptr[IDX(z  ,x,y  ,dim)] +
                             ptr[IDX(z  ,x,y+1,dim)] +
                             ptr[IDX(z+1,x,y  ,dim)] +
                             ptr[IDX(z+1,x,y+1,dim)])) );
};

template <const int BDIMX>
__device__ inline
float cell_coeff_BL_shfl (const float* __restrict__ ptr_gmem,
                          const int   z, 
                          const int   x, 
                          const int   y, 
                          const dim_t dim)
{
    float current = ptr_gmem[IDX(z,x,y,dim)];
    float right1  = __shfl_down(current, 1);

    /* For threads without neighbors: */
    if (threadIdx.x >= BDIMX-1 /* 1 */) right1 = ptr_gmem[IDX(z+1,x,y,dim)];

    float current_front = ptr_gmem[IDX(z,x,y+1,dim)];
    float right1_front  = __shfl_down(current_front, 1);

    /* For threads without neighbors: */
    if (threadIdx.x >= BDIMX-1 /* 1 */) right1_front = ptr_gmem[IDX(z+1,x,y+1,dim)];

    return ( 1.0f / ( 2.5f *(current       +
                             current_front +
                             right1        +
                             right1_front)) );
};

__device__ inline
float cell_coeff_TR (const float* __restrict__ ptr, 
                     const int   z, 
                     const int   x, 
                     const int   y, 
                     const dim_t dim)
{
    return ( 1.0f / ( 2.5f *(ptr[IDX(z  , x  , y  ,dim)] +
                             ptr[IDX(z  , x+1, y  ,dim)] +
                             ptr[IDX(z  , x  , y+1,dim)] +
                             ptr[IDX(z  , x+1, y+1,dim)])));
};

template <const int NX,
          const int NY,
          const int BDIMY>
__device__ inline
float cell_coeff_TR_smem (      float ptr_smem[NY][NX],
                          const float* __restrict__ ptr_gmem, 
                          const int   z, 
                          const int   x, 
                          const int   y, 
                          const dim_t dim)
{
    const int tx = threadIdx.x;
    const int ty = threadIdx.y;
    ///////////// intra-block communication///////////////////
    ptr_smem[ty][tx] = (z<dim.zsize && x<dim.xsize) ? ptr_gmem[IDX(z,x,y,dim)] : 0.0f;
    if (ty < 1) ptr_smem[ty+BDIMY][tx] = (z+dim.zsize && x+BDIMY<dim.xsize) ? ptr_gmem[IDX(z,x+BDIMY,y,dim)] : 0.0f;

    __syncthreads();

    const float current      = ptr_smem[ty  ][tx];
    const float current_down = ptr_smem[ty+1][tx];

    __syncthreads();

    ptr_smem[ty][tx] = (z<dim.zsize && x<dim.xsize) ? ptr_gmem[IDX(z,x,y+1,dim)] : 0.0f;
    if (ty < 1) ptr_smem[ty+BDIMY][tx] = (z+dim.zsize && x+BDIMY<dim.xsize) ? ptr_gmem[IDX(z,x+BDIMY,y+1,dim)] : 0.0f;

    __syncthreads();

    const float front      = ptr_smem[ty  ][tx];
    const float front_down = ptr_smem[ty+1][tx];

    return ( 1.0f / ( 2.5f *(current      +
                             current_down +
                             front        +
                             front_down)) );
};

__device__ inline
float cell_coeff_ARTM_BR(const float* __restrict__ ptr, 
                         const int   z, 
                         const int   x, 
                         const int   y, 
                         const dim_t dim)
{
    return ((1.0f / ptr[IDX(z  ,x  ,y,dim)]  +
             1.0f / ptr[IDX(z  ,x+1,y,dim)]  +
             1.0f / ptr[IDX(z+1,x  ,y,dim)]  +
             1.0f / ptr[IDX(z+1,x+1,y,dim)]) * 0.25f);
};

template <const int NX,
          const int NY,
          const int BDIMX,
          const int BDIMY>
__device__ inline
float cell_coeff_ARTM_BR_smem(      float ptr_smem[NY][NX],
                              const float* __restrict__ ptr_gmem, 
                              const int   z, 
                              const int   x, 
                              const int   y, 
                              const dim_t dim)
{
    const int tx = threadIdx.x;
    const int ty = threadIdx.y;
    ///////////// intra-block communication///////////////////
    ptr_smem[ty][tx] = ptr_gmem[IDX(z,x,y,dim)];
    if (tx < 1) ptr_smem[ty][tx+BDIMX] = ptr_gmem[IDX(z+BDIMX,x,y,dim)];
    if (ty < 1) ptr_smem[ty+BDIMY][tx] = ptr_gmem[IDX(z,x+BDIMY,y,dim)];
    if (tx == 0 && ty == 0) 
                ptr_smem[ty+BDIMY][tx+BDIMX] = ptr_gmem[IDX(z+BDIMX,x+BDIMY,y,dim)];
    __syncthreads();

    return ((1.0f / ptr_smem[ty  ][tx  ]  +
             1.0f / ptr_smem[ty+1][tx  ]  +
             1.0f / ptr_smem[ty  ][tx+1]  +
             1.0f / ptr_smem[ty+1][tx+1]) * 0.25f);
};

__device__ inline
float cell_coeff_ARTM_TL( const float* __restrict__ ptr, 
                         const int   z, 
                         const int   x, 
                         const int   y, 
                         const dim_t dim)
{
    return (1.0f / ptr[IDX(z,x,y,dim)]);
};

__device__ inline
float cell_coeff_ARTM_BL(const float* __restrict__ ptr, 
                         const int   z, 
                         const int   x, 
                         const int   y, 
                         const dim_t dim)
{
    return ((1.0f / ptr[IDX(z  ,x,y  ,dim)]  +
             1.0f / ptr[IDX(z  ,x,y+1,dim)]  +
             1.0f / ptr[IDX(z+1,x,y  ,dim)]  +
             1.0f / ptr[IDX(z+1,x,y+1,dim)]) * 0.25f);
};

template <const int BDIMX>
__device__ inline
float cell_coeff_ARTM_BL_shfl(const float* __restrict__ ptr_gmem, 
                              const int   z, 
                              const int   x, 
                              const int   y, 
                              const dim_t dim)
{
    float current = ((z < dim.zsize && x < dim.xsize) ? ptr_gmem[IDX(z,x,y,dim)] : 1.0f);
    float right1  = __shfl_down(current, 1);

    /* For threads without neighbors: */
    if (threadIdx.x >= BDIMX-1 /* 1 */) right1 = ((z+1 < dim.zsize && x < dim.xsize) ? ptr_gmem[IDX(z+1,x,y,dim)] : 1.0f);

    float current_front = ((z < dim.zsize && x < dim.xsize) ? ptr_gmem[IDX(z,x,y+1,dim)] : 1.0f);
    float right1_front  = __shfl_down(current_front, 1);

    /* For threads without neighbors: */
    if (threadIdx.x >= BDIMX-1 /* 1 */) right1_front = ((z+1 < dim.zsize && x < dim.xsize) ? ptr_gmem[IDX(z+1,x,y+1,dim)] : 1.0f);


    return ((1.0f / current       +
             1.0f / current_front +
             1.0f / right1        +
             1.0f / right1_front   ) * 0.25f);
};

__device__ inline
float cell_coeff_ARTM_TR(const float* __restrict__ ptr, 
                         const int   z, 
                         const int   x, 
                         const int   y, 
                         const dim_t dim)
{
    return ((1.0f / ptr[IDX(z,x  ,y  ,dim)]  +
             1.0f / ptr[IDX(z,x+1,y  ,dim)]  +
             1.0f / ptr[IDX(z,x  ,y+1,dim)]  +
             1.0f / ptr[IDX(z,x+1,y+1,dim)]) * 0.25f);
};

template <const int NX,
          const int NY,
          const int BDIMY>
__device__ inline
float cell_coeff_ARTM_TR_smem(      float ptr_smem[NY][NX],
                              const float* __restrict__ ptr_gmem, 
                              const int   z, 
                              const int   x, 
                              const int   y, 
                              const dim_t dim)
{
    const int tx = threadIdx.x;
    const int ty = threadIdx.y;
    ///////////// intra-block communication///////////////////
    ptr_smem[ty][tx] = ptr_gmem[IDX(z,x,y,dim)];
    if (ty < 1) ptr_smem[ty+BDIMY][tx] = ptr_gmem[IDX(z,x+BDIMY,y,dim)];
    __syncthreads();

    const float current      = ptr_smem[ty  ][tx];
    const float current_down = ptr_smem[ty+1][tx];

    __syncthreads();
    ptr_smem[ty][tx] = ptr_gmem[IDX(z,x,y+1,dim)];
    if (ty < 1) ptr_smem[ty+BDIMY][tx] = ptr_gmem[IDX(z,x+BDIMY,y+1,dim)];
    __syncthreads();

    const float front      = ptr_smem[ty  ][tx];
    const float front_down = ptr_smem[ty+1][tx];

    return ((1.0f / current      +
             1.0f / current_down +
             1.0f / front        +
             1.0f / front_down ) * 0.25f);
};

#ifdef SCELL_TR
template <const int HALO  = 4, 
          const int BDIMX = 32,
          const int BDIMY =  4,
          const int NX    = BDIMX+2*HALO,
          const int NY    = BDIMY+2*HALO>
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
                                   const dim_t  dim)
{
    const int z = blockIdx.x * blockDim.x + threadIdx.x + nz0; 
    const int x = blockIdx.y * blockDim.y + threadIdx.y + nx0;       

    __shared__ float bsmem[NY][NX];

    for(int y = ny0; y < nyf; y++)
    {
        const float c11 = cell_coeff_TR_smem      <NX,NY,BDIMY>(bsmem, cc11, z, x, y, dim);
        const float c12 = cell_coeff_TR_smem      <NX,NY,BDIMY>(bsmem, cc12, z, x, y, dim);
        const float c13 = cell_coeff_TR_smem      <NX,NY,BDIMY>(bsmem, cc13, z, x, y, dim);
        const float c14 = cell_coeff_ARTM_TR_smem <NX,NY,BDIMY>(bsmem, cc14, z, x, y, dim);
        const float c15 = cell_coeff_ARTM_TR_smem <NX,NY,BDIMY>(bsmem, cc15, z, x, y, dim);
        const float c16 = cell_coeff_ARTM_TR_smem <NX,NY,BDIMY>(bsmem, cc16, z, x, y, dim);
        const float c22 = cell_coeff_TR_smem      <NX,NY,BDIMY>(bsmem, cc22, z, x, y, dim);
        const float c23 = cell_coeff_TR_smem      <NX,NY,BDIMY>(bsmem, cc23, z, x, y, dim);
        const float c24 = cell_coeff_ARTM_TR_smem <NX,NY,BDIMY>(bsmem, cc24, z, x, y, dim);
        const float c25 = cell_coeff_ARTM_TR_smem <NX,NY,BDIMY>(bsmem, cc25, z, x, y, dim);
        const float c26 = cell_coeff_ARTM_TR_smem <NX,NY,BDIMY>(bsmem, cc26, z, x, y, dim);
        const float c33 = cell_coeff_TR_smem      <NX,NY,BDIMY>(bsmem, cc33, z, x, y, dim);
        const float c34 = cell_coeff_ARTM_TR_smem <NX,NY,BDIMY>(bsmem, cc34, z, x, y, dim);
        const float c35 = cell_coeff_ARTM_TR_smem <NX,NY,BDIMY>(bsmem, cc35, z, x, y, dim);
        const float c36 = cell_coeff_ARTM_TR_smem <NX,NY,BDIMY>(bsmem, cc36, z, x, y, dim);
        const float c44 = cell_coeff_TR_smem      <NX,NY,BDIMY>(bsmem, cc44, z, x, y, dim);
        const float c45 = cell_coeff_ARTM_TR_smem <NX,NY,BDIMY>(bsmem, cc45, z, x, y, dim);
        const float c46 = cell_coeff_ARTM_TR_smem <NX,NY,BDIMY>(bsmem, cc46, z, x, y, dim);
        const float c55 = cell_coeff_TR_smem      <NX,NY,BDIMY>(bsmem, cc55, z, x, y, dim);
        const float c56 = cell_coeff_ARTM_TR_smem <NX,NY,BDIMY>(bsmem, cc56, z, x, y, dim);
        const float c66 = cell_coeff_TR_smem      <NX,NY,BDIMY>(bsmem, cc66, z, x, y, dim);

        const float u_x = stencil_X_smem <NX,NY,HALO,BDIMY>(SX, bsmem, vxu, dxi, z, x, y, dim);
        const float v_x = stencil_X_smem <NX,NY,HALO,BDIMY>(SX, bsmem, vxv, dxi, z, x, y, dim);
        const float w_x = stencil_X_smem <NX,NY,HALO,BDIMY>(SX, bsmem, vxw, dxi, z, x, y, dim);
        
        float u_y, v_y, w_y;
        if (z < nzf && x < nxf)
        {
            u_y = stencil_Y (SY, vyu, dyi, z, x, y, dim);
            v_y = stencil_Y (SY, vyv, dyi, z, x, y, dim);
            w_y = stencil_Y (SY, vyw, dyi, z, x, y, dim);
        }
        
        const float u_z = stencil_Z_shfl <HALO,BDIMX>(SZ, vzu, dzi, z, x, y, dim);
        const float v_z = stencil_Z_shfl <HALO,BDIMX>(SZ, vzv, dzi, z, x, y, dim);
        const float w_z = stencil_Z_shfl <HALO,BDIMX>(SZ, vzw, dzi, z, x, y, dim);
        
        if (z < nzf && x < nxf)
        {
            stress_update (sxxptr,c11,c12,c13,c14,c15,c16,z,x,y,dt,u_x,u_y,u_z,v_x,v_y,v_z,w_x,w_y,w_z,dim );
            stress_update (syyptr,c12,c22,c23,c24,c25,c26,z,x,y,dt,u_x,u_y,u_z,v_x,v_y,v_z,w_x,w_y,w_z,dim );
            stress_update (szzptr,c13,c23,c33,c34,c35,c36,z,x,y,dt,u_x,u_y,u_z,v_x,v_y,v_z,w_x,w_y,w_z,dim );
            stress_update (syzptr,c14,c24,c34,c44,c45,c46,z,x,y,dt,u_x,u_y,u_z,v_x,v_y,v_z,w_x,w_y,w_z,dim );
            stress_update (sxzptr,c15,c25,c35,c45,c55,c56,z,x,y,dt,u_x,u_y,u_z,v_x,v_y,v_z,w_x,w_y,w_z,dim );
            stress_update (sxyptr,c16,c26,c36,c46,c56,c66,z,x,y,dt,u_x,u_y,u_z,v_x,v_y,v_z,w_x,w_y,w_z,dim );
        }
    }
}
#else
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
                                   const dim_t  dim)
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
                const float c11 = cell_coeff_TR      (cc11, z, x, y, dim);
                const float c12 = cell_coeff_TR      (cc12, z, x, y, dim);
                const float c13 = cell_coeff_TR      (cc13, z, x, y, dim);
                const float c14 = cell_coeff_ARTM_TR (cc14, z, x, y, dim);
                const float c15 = cell_coeff_ARTM_TR (cc15, z, x, y, dim);
                const float c16 = cell_coeff_ARTM_TR (cc16, z, x, y, dim);
                const float c22 = cell_coeff_TR      (cc22, z, x, y, dim);
                const float c23 = cell_coeff_TR      (cc23, z, x, y, dim);
                const float c24 = cell_coeff_ARTM_TR (cc24, z, x, y, dim);
                const float c25 = cell_coeff_ARTM_TR (cc25, z, x, y, dim);
                const float c26 = cell_coeff_ARTM_TR (cc26, z, x, y, dim);
                const float c33 = cell_coeff_TR      (cc33, z, x, y, dim);
                const float c34 = cell_coeff_ARTM_TR (cc34, z, x, y, dim);
                const float c35 = cell_coeff_ARTM_TR (cc35, z, x, y, dim);
                const float c36 = cell_coeff_ARTM_TR (cc36, z, x, y, dim);
                const float c44 = cell_coeff_TR      (cc44, z, x, y, dim);
                const float c45 = cell_coeff_ARTM_TR (cc45, z, x, y, dim);
                const float c46 = cell_coeff_ARTM_TR (cc46, z, x, y, dim);
                const float c55 = cell_coeff_TR      (cc55, z, x, y, dim);
                const float c56 = cell_coeff_ARTM_TR (cc56, z, x, y, dim);
                const float c66 = cell_coeff_TR      (cc66, z, x, y, dim);
                
                const float u_x = stencil_X (SX, vxu, dxi, z, x, y, dim);
                const float v_x = stencil_X (SX, vxv, dxi, z, x, y, dim);
                const float w_x = stencil_X (SX, vxw, dxi, z, x, y, dim);
                
                const float u_y = stencil_Y (SY, vyu, dyi, z, x, y, dim);
                const float v_y = stencil_Y (SY, vyv, dyi, z, x, y, dim);
                const float w_y = stencil_Y (SY, vyw, dyi, z, x, y, dim);
                
                const float u_z = stencil_Z (SZ, vzu, dzi, z, x, y, dim);
                const float v_z = stencil_Z (SZ, vzv, dzi, z, x, y, dim);
                const float w_z = stencil_Z (SZ, vzw, dzi, z, x, y, dim);
                
                stress_update (sxxptr,c11,c12,c13,c14,c15,c16,z,x,y,dt,u_x,u_y,u_z,v_x,v_y,v_z,w_x,w_y,w_z,dim );
                stress_update (syyptr,c12,c22,c23,c24,c25,c26,z,x,y,dt,u_x,u_y,u_z,v_x,v_y,v_z,w_x,w_y,w_z,dim );
                stress_update (szzptr,c13,c23,c33,c34,c35,c36,z,x,y,dt,u_x,u_y,u_z,v_x,v_y,v_z,w_x,w_y,w_z,dim );
                stress_update (syzptr,c14,c24,c34,c44,c45,c46,z,x,y,dt,u_x,u_y,u_z,v_x,v_y,v_z,w_x,w_y,w_z,dim );
                stress_update (sxzptr,c15,c25,c35,c45,c55,c56,z,x,y,dt,u_x,u_y,u_z,v_x,v_y,v_z,w_x,w_y,w_z,dim );
                stress_update (sxyptr,c16,c26,c36,c46,c56,c66,z,x,y,dt,u_x,u_y,u_z,v_x,v_y,v_z,w_x,w_y,w_z,dim );
            }
        }
    }
}
#endif

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
                                 const dim_t  dim,
                                 void*        stream)
{
    const int block_dim_x = 32;
    const int block_dim_y = 4;

    dim3 grid_dim( ((nzf-nz0) + block_dim_x-1)/block_dim_x,
                   ((nxf-nx0) + block_dim_y-1)/block_dim_y,
                      1 );
    dim3 block_dim(block_dim_x, block_dim_y, 1);

    cudaStream_t s = (cudaStream_t) stream;

#ifdef SCELL_TR
    const int HALO = 4;
    compute_component_scell_TR_cuda_k<HALO,block_dim_x,block_dim_y><<<grid_dim, block_dim, 0, s>>>
        (sxxptr, syyptr, szzptr, syzptr, sxzptr, sxyptr,
         vxu, vxv, vxw, vyu, vyv, vyw, vzu, vzv, vzw,
         cc11, cc12, cc13, cc14, cc15, cc16,
         cc22, cc23, cc24, cc25, cc26,
         cc33, cc34, cc35, cc36,
         cc44, cc45, cc46,
         cc55, cc56,
         cc66,
         dt, dzi, dxi, dyi, nz0, nzf, nx0, nxf, ny0, nyf, SZ, SX, SY, dim);
#else
    compute_component_scell_TR_cuda_k<<<grid_dim, block_dim, 0, s>>>
        (sxxptr, syyptr, szzptr, syzptr, sxzptr, sxyptr,
         vxu, vxv, vxw, vyu, vyv, vyw, vzu, vzv, vzw,
         cc11, cc12, cc13, cc14, cc15, cc16,
         cc22, cc23, cc24, cc25, cc26,
         cc33, cc34, cc35, cc36,
         cc44, cc45, cc46,
         cc55, cc56,
         cc66,
         dt, dzi, dxi, dyi, nz0, nzf, nx0, nxf, ny0, nyf, SZ, SX, SY, dim);
#endif 
    CUDA_CHECK(cudaGetLastError());
};

#ifdef SCELL_TL
template <const int HALO  = 4, 
          const int BDIMX = 32,
          const int BDIMY =  4,
          const int NX    = BDIMX+2*HALO,
          const int NY    = BDIMY+2*HALO>
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
                                   const dim_t  dim)
{
    const int z = blockIdx.x * blockDim.x + threadIdx.x + nz0; 
    const int x = blockIdx.y * blockDim.y + threadIdx.y + nx0;       

    __shared__ float bsmem[NY][NX];

    for(int y = ny0; y < nyf; y++)
    {
        float c11, c12, c13, c14, c15, c16;
        float c22, c23, c24, c25, c26;
        float c33, c34, c35, c36;
        float c44, c45, c46;
        float c55, c56;
        float c66;

        if (z < dim.zsize && x < dim.xsize)
        {
            c11 = cell_coeff_TL      (cc11, z, x, y, dim);
            c12 = cell_coeff_TL      (cc12, z, x, y, dim);
            c13 = cell_coeff_TL      (cc13, z, x, y, dim);
            c14 = cell_coeff_ARTM_TL (cc14, z, x, y, dim);
            c15 = cell_coeff_ARTM_TL (cc15, z, x, y, dim);
            c16 = cell_coeff_ARTM_TL (cc16, z, x, y, dim);
            c22 = cell_coeff_TL      (cc22, z, x, y, dim);
            c23 = cell_coeff_TL      (cc23, z, x, y, dim);
            c24 = cell_coeff_ARTM_TL (cc24, z, x, y, dim);
            c25 = cell_coeff_ARTM_TL (cc25, z, x, y, dim);
            c26 = cell_coeff_ARTM_TL (cc26, z, x, y, dim);
            c33 = cell_coeff_TL      (cc33, z, x, y, dim);
            c34 = cell_coeff_ARTM_TL (cc34, z, x, y, dim);
            c35 = cell_coeff_ARTM_TL (cc35, z, x, y, dim);
            c36 = cell_coeff_ARTM_TL (cc36, z, x, y, dim);
            c44 = cell_coeff_TL      (cc44, z, x, y, dim);
            c45 = cell_coeff_ARTM_TL (cc45, z, x, y, dim);
            c46 = cell_coeff_ARTM_TL (cc46, z, x, y, dim);
            c55 = cell_coeff_TL      (cc55, z, x, y, dim);
            c56 = cell_coeff_ARTM_TL (cc56, z, x, y, dim);
            c66 = cell_coeff_TL      (cc66, z, x, y, dim);
        }

        const float u_x = stencil_X_smem <NX,NY,HALO,BDIMY>(SX, bsmem, vxu, dxi, z, x, y, dim);
        const float v_x = stencil_X_smem <NX,NY,HALO,BDIMY>(SX, bsmem, vxv, dxi, z, x, y, dim);
        const float w_x = stencil_X_smem <NX,NY,HALO,BDIMY>(SX, bsmem, vxw, dxi, z, x, y, dim);
        
        float u_y, v_y, w_y;
        if (z < nzf && x < nxf)
        {
            u_y = stencil_Y (SY, vyu, dyi, z, x, y, dim);
            v_y = stencil_Y (SY, vyv, dyi, z, x, y, dim);
            w_y = stencil_Y (SY, vyw, dyi, z, x, y, dim);
        }
        
        const float u_z = stencil_Z_shfl <HALO,BDIMX>(SZ, vzu, dzi, z, x, y, dim);
        const float v_z = stencil_Z_shfl <HALO,BDIMX>(SZ, vzv, dzi, z, x, y, dim);
        const float w_z = stencil_Z_shfl <HALO,BDIMX>(SZ, vzw, dzi, z, x, y, dim);
        
        if (z < nzf && x < nxf)
        {
            stress_update (sxxptr,c11,c12,c13,c14,c15,c16,z,x,y,dt,u_x,u_y,u_z,v_x,v_y,v_z,w_x,w_y,w_z,dim );
            stress_update (syyptr,c12,c22,c23,c24,c25,c26,z,x,y,dt,u_x,u_y,u_z,v_x,v_y,v_z,w_x,w_y,w_z,dim );
            stress_update (szzptr,c13,c23,c33,c34,c35,c36,z,x,y,dt,u_x,u_y,u_z,v_x,v_y,v_z,w_x,w_y,w_z,dim );
            stress_update (syzptr,c14,c24,c34,c44,c45,c46,z,x,y,dt,u_x,u_y,u_z,v_x,v_y,v_z,w_x,w_y,w_z,dim );
            stress_update (sxzptr,c15,c25,c35,c45,c55,c56,z,x,y,dt,u_x,u_y,u_z,v_x,v_y,v_z,w_x,w_y,w_z,dim );
            stress_update (sxyptr,c16,c26,c36,c46,c56,c66,z,x,y,dt,u_x,u_y,u_z,v_x,v_y,v_z,w_x,w_y,w_z,dim );
        }
    }
}
#else
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
                                   const dim_t  dim)
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
                const float c11 = cell_coeff_TL      (cc11, z, x, y, dim);
                const float c12 = cell_coeff_TL      (cc12, z, x, y, dim);
                const float c13 = cell_coeff_TL      (cc13, z, x, y, dim);
                const float c14 = cell_coeff_ARTM_TL (cc14, z, x, y, dim);
                const float c15 = cell_coeff_ARTM_TL (cc15, z, x, y, dim);
                const float c16 = cell_coeff_ARTM_TL (cc16, z, x, y, dim);
                const float c22 = cell_coeff_TL      (cc22, z, x, y, dim);
                const float c23 = cell_coeff_TL      (cc23, z, x, y, dim);
                const float c24 = cell_coeff_ARTM_TL (cc24, z, x, y, dim);
                const float c25 = cell_coeff_ARTM_TL (cc25, z, x, y, dim);
                const float c26 = cell_coeff_ARTM_TL (cc26, z, x, y, dim);
                const float c33 = cell_coeff_TL      (cc33, z, x, y, dim);
                const float c34 = cell_coeff_ARTM_TL (cc34, z, x, y, dim);
                const float c35 = cell_coeff_ARTM_TL (cc35, z, x, y, dim);
                const float c36 = cell_coeff_ARTM_TL (cc36, z, x, y, dim);
                const float c44 = cell_coeff_TL      (cc44, z, x, y, dim);
                const float c45 = cell_coeff_ARTM_TL (cc45, z, x, y, dim);
                const float c46 = cell_coeff_ARTM_TL (cc46, z, x, y, dim);
                const float c55 = cell_coeff_TL      (cc55, z, x, y, dim);
                const float c56 = cell_coeff_ARTM_TL (cc56, z, x, y, dim);
                const float c66 = cell_coeff_TL      (cc66, z, x, y, dim);
                
                const float u_x = stencil_X (SX, vxu, dxi, z, x, y, dim);
                const float v_x = stencil_X (SX, vxv, dxi, z, x, y, dim);
                const float w_x = stencil_X (SX, vxw, dxi, z, x, y, dim);
                
                const float u_y = stencil_Y (SY, vyu, dyi, z, x, y, dim);
                const float v_y = stencil_Y (SY, vyv, dyi, z, x, y, dim);
                const float w_y = stencil_Y (SY, vyw, dyi, z, x, y, dim);
                
                const float u_z = stencil_Z (SZ, vzu, dzi, z, x, y, dim);
                const float v_z = stencil_Z (SZ, vzv, dzi, z, x, y, dim);
                const float w_z = stencil_Z (SZ, vzw, dzi, z, x, y, dim);
                
                stress_update (sxxptr,c11,c12,c13,c14,c15,c16,z,x,y,dt,u_x,u_y,u_z,v_x,v_y,v_z,w_x,w_y,w_z,dim );
                stress_update (syyptr,c12,c22,c23,c24,c25,c26,z,x,y,dt,u_x,u_y,u_z,v_x,v_y,v_z,w_x,w_y,w_z,dim );
                stress_update (szzptr,c13,c23,c33,c34,c35,c36,z,x,y,dt,u_x,u_y,u_z,v_x,v_y,v_z,w_x,w_y,w_z,dim );
                stress_update (syzptr,c14,c24,c34,c44,c45,c46,z,x,y,dt,u_x,u_y,u_z,v_x,v_y,v_z,w_x,w_y,w_z,dim );
                stress_update (sxzptr,c15,c25,c35,c45,c55,c56,z,x,y,dt,u_x,u_y,u_z,v_x,v_y,v_z,w_x,w_y,w_z,dim );
                stress_update (sxyptr,c16,c26,c36,c46,c56,c66,z,x,y,dt,u_x,u_y,u_z,v_x,v_y,v_z,w_x,w_y,w_z,dim );
            }
        }
    }
}
#endif

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
                                 const dim_t  dim,
                                 void*        stream)
{
    const int block_dim_x = 32;
    const int block_dim_y = 4;

    dim3 grid_dim( ((nzf-nz0) + block_dim_x-1)/block_dim_x,
                   ((nxf-nx0) + block_dim_y-1)/block_dim_y,
                      1 );
    dim3 block_dim(block_dim_x, block_dim_y, 1);

    cudaStream_t s = (cudaStream_t) stream;

#ifdef SCELL_TL
    const int HALO = 4;
    compute_component_scell_TL_cuda_k<HALO,block_dim_x,block_dim_y><<<grid_dim, block_dim, 0, s>>>
        (sxxptr, syyptr, szzptr, syzptr, sxzptr, sxyptr,
         vxu, vxv, vxw, vyu, vyv, vyw, vzu, vzv, vzw,
         cc11, cc12, cc13, cc14, cc15, cc16,
         cc22, cc23, cc24, cc25, cc26,
         cc33, cc34, cc35, cc36,
         cc44, cc45, cc46,
         cc55, cc56,
         cc66,
         dt, dzi, dxi, dyi, nz0, nzf, nx0, nxf, ny0, nyf, SZ, SX, SY, dim);
#else
    compute_component_scell_TL_cuda_k<<<grid_dim, block_dim, 0, s>>>
        (sxxptr, syyptr, szzptr, syzptr, sxzptr, sxyptr,
         vxu, vxv, vxw, vyu, vyv, vyw, vzu, vzv, vzw,
         cc11, cc12, cc13, cc14, cc15, cc16,
         cc22, cc23, cc24, cc25, cc26,
         cc33, cc34, cc35, cc36,
         cc44, cc45, cc46,
         cc55, cc56,
         cc66,
         dt, dzi, dxi, dyi, nz0, nzf, nx0, nxf, ny0, nyf, SZ, SX, SY, dim);
#endif 
    CUDA_CHECK(cudaGetLastError());
};

#ifdef SCELL_BR
template <const int HALO  = 4, 
          const int BDIMX = 32,
          const int BDIMY =  4,
          const int NX    = BDIMX+2*HALO,
          const int NY    = BDIMY+2*HALO>
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
                                   const dim_t  dim)
{
    const int z = blockIdx.x * blockDim.x + threadIdx.x + nz0; 
    const int x = blockIdx.y * blockDim.y + threadIdx.y + nx0;       

    __shared__ float bsmem[NY][NX];

    for(int y = ny0; y < nyf; y++)
    {
        const float c11 = cell_coeff_BR_smem      <NX,NY,BDIMX,BDIMY>(bsmem, cc11, z, x, y, dim);
        const float c12 = cell_coeff_BR_smem      <NX,NY,BDIMX,BDIMY>(bsmem, cc12, z, x, y, dim);
        const float c13 = cell_coeff_BR_smem      <NX,NY,BDIMX,BDIMY>(bsmem, cc13, z, x, y, dim);
        const float c14 = cell_coeff_ARTM_BR_smem <NX,NY,BDIMX,BDIMY>(bsmem, cc14, z, x, y, dim);
        const float c15 = cell_coeff_ARTM_BR_smem <NX,NY,BDIMX,BDIMY>(bsmem, cc15, z, x, y, dim);
        const float c16 = cell_coeff_ARTM_BR_smem <NX,NY,BDIMX,BDIMY>(bsmem, cc16, z, x, y, dim);
        const float c22 = cell_coeff_BR_smem      <NX,NY,BDIMX,BDIMY>(bsmem, cc22, z, x, y, dim);
        const float c23 = cell_coeff_BR_smem      <NX,NY,BDIMX,BDIMY>(bsmem, cc23, z, x, y, dim);
        const float c24 = cell_coeff_ARTM_BR_smem <NX,NY,BDIMX,BDIMY>(bsmem, cc24, z, x, y, dim);
        const float c25 = cell_coeff_ARTM_BR_smem <NX,NY,BDIMX,BDIMY>(bsmem, cc25, z, x, y, dim);
        const float c26 = cell_coeff_ARTM_BR_smem <NX,NY,BDIMX,BDIMY>(bsmem, cc26, z, x, y, dim);
        const float c33 = cell_coeff_BR_smem      <NX,NY,BDIMX,BDIMY>(bsmem, cc33, z, x, y, dim);
        const float c34 = cell_coeff_ARTM_BR_smem <NX,NY,BDIMX,BDIMY>(bsmem, cc34, z, x, y, dim);
        const float c35 = cell_coeff_ARTM_BR_smem <NX,NY,BDIMX,BDIMY>(bsmem, cc35, z, x, y, dim);
        const float c36 = cell_coeff_ARTM_BR_smem <NX,NY,BDIMX,BDIMY>(bsmem, cc36, z, x, y, dim);
        const float c44 = cell_coeff_BR_smem      <NX,NY,BDIMX,BDIMY>(bsmem, cc44, z, x, y, dim);
        const float c45 = cell_coeff_ARTM_BR_smem <NX,NY,BDIMX,BDIMY>(bsmem, cc45, z, x, y, dim);
        const float c46 = cell_coeff_ARTM_BR_smem <NX,NY,BDIMX,BDIMY>(bsmem, cc46, z, x, y, dim);
        const float c55 = cell_coeff_BR_smem      <NX,NY,BDIMX,BDIMY>(bsmem, cc55, z, x, y, dim);
        const float c56 = cell_coeff_ARTM_BR_smem <NX,NY,BDIMX,BDIMY>(bsmem, cc56, z, x, y, dim);
        const float c66 = cell_coeff_BR_smem      <NX,NY,BDIMX,BDIMY>(bsmem, cc66, z, x, y, dim);

        const float u_x = stencil_X_smem <NX,NY,HALO,BDIMY>(SX, bsmem, vxu, dxi, z, x, y, dim);
        const float v_x = stencil_X_smem <NX,NY,HALO,BDIMY>(SX, bsmem, vxv, dxi, z, x, y, dim);
        const float w_x = stencil_X_smem <NX,NY,HALO,BDIMY>(SX, bsmem, vxw, dxi, z, x, y, dim);
        
        float u_y, v_y, w_y;
        if (z < nzf && x < nxf)
        {
            u_y = stencil_Y (SY, vyu, dyi, z, x, y, dim);
            v_y = stencil_Y (SY, vyv, dyi, z, x, y, dim);
            w_y = stencil_Y (SY, vyw, dyi, z, x, y, dim);
        }

        const float u_z = stencil_Z_shfl <HALO,BDIMX>(SZ, vzu, dzi, z, x, y, dim);
        const float v_z = stencil_Z_shfl <HALO,BDIMX>(SZ, vzv, dzi, z, x, y, dim);
        const float w_z = stencil_Z_shfl <HALO,BDIMX>(SZ, vzw, dzi, z, x, y, dim);
        
        if (z < nzf && x < nxf)
        {
            stress_update (sxxptr,c11,c12,c13,c14,c15,c16,z,x,y,dt,u_x,u_y,u_z,v_x,v_y,v_z,w_x,w_y,w_z,dim );
            stress_update (syyptr,c12,c22,c23,c24,c25,c26,z,x,y,dt,u_x,u_y,u_z,v_x,v_y,v_z,w_x,w_y,w_z,dim );
            stress_update (szzptr,c13,c23,c33,c34,c35,c36,z,x,y,dt,u_x,u_y,u_z,v_x,v_y,v_z,w_x,w_y,w_z,dim );
            stress_update (syzptr,c14,c24,c34,c44,c45,c46,z,x,y,dt,u_x,u_y,u_z,v_x,v_y,v_z,w_x,w_y,w_z,dim );
            stress_update (sxzptr,c15,c25,c35,c45,c55,c56,z,x,y,dt,u_x,u_y,u_z,v_x,v_y,v_z,w_x,w_y,w_z,dim );
            stress_update (sxyptr,c16,c26,c36,c46,c56,c66,z,x,y,dt,u_x,u_y,u_z,v_x,v_y,v_z,w_x,w_y,w_z,dim );
        }
    }
}
#else
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
                                   const dim_t  dim)
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
                const float c11 = cell_coeff_BR      (cc11, z, x, y, dim);
                const float c12 = cell_coeff_BR      (cc12, z, x, y, dim);
                const float c13 = cell_coeff_BR      (cc13, z, x, y, dim);
                const float c14 = cell_coeff_ARTM_BR (cc14, z, x, y, dim);
                const float c15 = cell_coeff_ARTM_BR (cc15, z, x, y, dim);
                const float c16 = cell_coeff_ARTM_BR (cc16, z, x, y, dim);
                const float c22 = cell_coeff_BR      (cc22, z, x, y, dim);
                const float c23 = cell_coeff_BR      (cc23, z, x, y, dim);
                const float c24 = cell_coeff_ARTM_BR (cc24, z, x, y, dim);
                const float c25 = cell_coeff_ARTM_BR (cc25, z, x, y, dim);
                const float c26 = cell_coeff_ARTM_BR (cc26, z, x, y, dim);
                const float c33 = cell_coeff_BR      (cc33, z, x, y, dim);
                const float c34 = cell_coeff_ARTM_BR (cc34, z, x, y, dim);
                const float c35 = cell_coeff_ARTM_BR (cc35, z, x, y, dim);
                const float c36 = cell_coeff_ARTM_BR (cc36, z, x, y, dim);
                const float c44 = cell_coeff_BR      (cc44, z, x, y, dim);
                const float c45 = cell_coeff_ARTM_BR (cc45, z, x, y, dim);
                const float c46 = cell_coeff_ARTM_BR (cc46, z, x, y, dim);
                const float c55 = cell_coeff_BR      (cc55, z, x, y, dim);
                const float c56 = cell_coeff_ARTM_BR (cc56, z, x, y, dim);
                const float c66 = cell_coeff_BR      (cc66, z, x, y, dim);
                
                const float u_x = stencil_X (SX, vxu, dxi, z, x, y, dim);
                const float v_x = stencil_X (SX, vxv, dxi, z, x, y, dim);
                const float w_x = stencil_X (SX, vxw, dxi, z, x, y, dim);
                
                const float u_y = stencil_Y (SY, vyu, dyi, z, x, y, dim);
                const float v_y = stencil_Y (SY, vyv, dyi, z, x, y, dim);
                const float w_y = stencil_Y (SY, vyw, dyi, z, x, y, dim);
                
                const float u_z = stencil_Z (SZ, vzu, dzi, z, x, y, dim);
                const float v_z = stencil_Z (SZ, vzv, dzi, z, x, y, dim);
                const float w_z = stencil_Z (SZ, vzw, dzi, z, x, y, dim);
                
                stress_update (sxxptr,c11,c12,c13,c14,c15,c16,z,x,y,dt,u_x,u_y,u_z,v_x,v_y,v_z,w_x,w_y,w_z,dim );
                stress_update (syyptr,c12,c22,c23,c24,c25,c26,z,x,y,dt,u_x,u_y,u_z,v_x,v_y,v_z,w_x,w_y,w_z,dim );
                stress_update (szzptr,c13,c23,c33,c34,c35,c36,z,x,y,dt,u_x,u_y,u_z,v_x,v_y,v_z,w_x,w_y,w_z,dim );
                stress_update (syzptr,c14,c24,c34,c44,c45,c46,z,x,y,dt,u_x,u_y,u_z,v_x,v_y,v_z,w_x,w_y,w_z,dim );
                stress_update (sxzptr,c15,c25,c35,c45,c55,c56,z,x,y,dt,u_x,u_y,u_z,v_x,v_y,v_z,w_x,w_y,w_z,dim );
                stress_update (sxyptr,c16,c26,c36,c46,c56,c66,z,x,y,dt,u_x,u_y,u_z,v_x,v_y,v_z,w_x,w_y,w_z,dim );
            }
        }
    }
}
#endif

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
                                 const dim_t  dim,
                                 void*        stream)
{
    const int block_dim_x = 32;
    const int block_dim_y = 4;

    dim3 grid_dim( ((nzf-nz0) + block_dim_x-1)/block_dim_x,
                   ((nxf-nx0) + block_dim_y-1)/block_dim_y,
                      1 );
    dim3 block_dim(block_dim_x, block_dim_y, 1);

    cudaStream_t s = (cudaStream_t) stream;

#ifdef SCELL_BR
    const int HALO = 4;
    compute_component_scell_BR_cuda_k<HALO,block_dim_x,block_dim_y><<<grid_dim, block_dim, 0, s>>>
        (sxxptr, syyptr, szzptr, syzptr, sxzptr, sxyptr,
         vxu, vxv, vxw, vyu, vyv, vyw, vzu, vzv, vzw,
         cc11, cc12, cc13, cc14, cc15, cc16,
         cc22, cc23, cc24, cc25, cc26,
         cc33, cc34, cc35, cc36,
         cc44, cc45, cc46,
         cc55, cc56,
         cc66,
         dt, dzi, dxi, dyi, nz0, nzf, nx0, nxf, ny0, nyf, SZ, SX, SY, dim);
#else
    compute_component_scell_BR_cuda_k<<<grid_dim, block_dim, 0, s>>>
        (sxxptr, syyptr, szzptr, syzptr, sxzptr, sxyptr,
         vxu, vxv, vxw, vyu, vyv, vyw, vzu, vzv, vzw,
         cc11, cc12, cc13, cc14, cc15, cc16,
         cc22, cc23, cc24, cc25, cc26,
         cc33, cc34, cc35, cc36,
         cc44, cc45, cc46,
         cc55, cc56,
         cc66,
         dt, dzi, dxi, dyi, nz0, nzf, nx0, nxf, ny0, nyf, SZ, SX, SY, dim);
#endif 
    CUDA_CHECK(cudaGetLastError());
};


#ifdef SCELL_BL
template <const int HALO  = 4, 
          const int BDIMX = 32,
          const int BDIMY =  4,
          const int NX    = BDIMX+2*HALO,
          const int NY    = BDIMY+2*HALO>
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
                                   const dim_t  dim)
{
    const int z = blockIdx.x * blockDim.x + threadIdx.x + nz0; 
    const int x = blockIdx.y * blockDim.y + threadIdx.y + nx0;       

    __shared__ float bsmem[NY][NX];

    for(int y = ny0; y < nyf; y++)
    {
        float c11, c12, c13, c14, c15, c16;
        float c22, c23, c24, c25, c26;
        float c33, c34, c35, c36;
        float c44, c45, c46;
        float c55, c56;
        float c66;

        if (z < dim.zsize && x < dim.xsize)
        {
            c11 = cell_coeff_BL_shfl      <BDIMX>(cc11, z, x, y, dim);
            c12 = cell_coeff_BL_shfl      <BDIMX>(cc12, z, x, y, dim);
            c13 = cell_coeff_BL_shfl      <BDIMX>(cc13, z, x, y, dim);
            c14 = cell_coeff_ARTM_BL_shfl <BDIMX>(cc14, z, x, y, dim);
            c15 = cell_coeff_ARTM_BL_shfl <BDIMX>(cc15, z, x, y, dim);
            c16 = cell_coeff_ARTM_BL_shfl <BDIMX>(cc16, z, x, y, dim);
            c22 = cell_coeff_BL_shfl      <BDIMX>(cc22, z, x, y, dim);
            c23 = cell_coeff_BL_shfl      <BDIMX>(cc23, z, x, y, dim);
            c24 = cell_coeff_ARTM_BL_shfl <BDIMX>(cc24, z, x, y, dim);
            c25 = cell_coeff_ARTM_BL_shfl <BDIMX>(cc25, z, x, y, dim);
            c26 = cell_coeff_ARTM_BL_shfl <BDIMX>(cc26, z, x, y, dim);
            c33 = cell_coeff_BL_shfl      <BDIMX>(cc33, z, x, y, dim);
            c34 = cell_coeff_ARTM_BL_shfl <BDIMX>(cc34, z, x, y, dim);
            c35 = cell_coeff_ARTM_BL_shfl <BDIMX>(cc35, z, x, y, dim);
            c36 = cell_coeff_ARTM_BL_shfl <BDIMX>(cc36, z, x, y, dim);
            c44 = cell_coeff_BL_shfl      <BDIMX>(cc44, z, x, y, dim);
            c45 = cell_coeff_ARTM_BL_shfl <BDIMX>(cc45, z, x, y, dim);
            c46 = cell_coeff_ARTM_BL_shfl <BDIMX>(cc46, z, x, y, dim);
            c55 = cell_coeff_BL_shfl      <BDIMX>(cc55, z, x, y, dim);
            c56 = cell_coeff_ARTM_BL_shfl <BDIMX>(cc56, z, x, y, dim);
            c66 = cell_coeff_BL_shfl      <BDIMX>(cc66, z, x, y, dim);
        }

        const float u_x = stencil_X_smem <NX,NY,HALO,BDIMY>(SX, bsmem, vxu, dxi, z, x, y, dim);
        const float v_x = stencil_X_smem <NX,NY,HALO,BDIMY>(SX, bsmem, vxv, dxi, z, x, y, dim);
        const float w_x = stencil_X_smem <NX,NY,HALO,BDIMY>(SX, bsmem, vxw, dxi, z, x, y, dim);
        
        float u_y, v_y, w_y;
        if (z < nzf && x < nxf)
        {
            u_y = stencil_Y (SY, vyu, dyi, z, x, y, dim);
            v_y = stencil_Y (SY, vyv, dyi, z, x, y, dim);
            w_y = stencil_Y (SY, vyw, dyi, z, x, y, dim);
        }
        
        const float u_z = stencil_Z_shfl <HALO,BDIMX>(SZ, vzu, dzi, z, x, y, dim);
        const float v_z = stencil_Z_shfl <HALO,BDIMX>(SZ, vzv, dzi, z, x, y, dim);
        const float w_z = stencil_Z_shfl <HALO,BDIMX>(SZ, vzw, dzi, z, x, y, dim);
        
        if (z < nzf && x < nxf)
        {
            stress_update (sxxptr,c11,c12,c13,c14,c15,c16,z,x,y,dt,u_x,u_y,u_z,v_x,v_y,v_z,w_x,w_y,w_z,dim );
            stress_update (syyptr,c12,c22,c23,c24,c25,c26,z,x,y,dt,u_x,u_y,u_z,v_x,v_y,v_z,w_x,w_y,w_z,dim );
            stress_update (szzptr,c13,c23,c33,c34,c35,c36,z,x,y,dt,u_x,u_y,u_z,v_x,v_y,v_z,w_x,w_y,w_z,dim );
            stress_update (syzptr,c14,c24,c34,c44,c45,c46,z,x,y,dt,u_x,u_y,u_z,v_x,v_y,v_z,w_x,w_y,w_z,dim );
            stress_update (sxzptr,c15,c25,c35,c45,c55,c56,z,x,y,dt,u_x,u_y,u_z,v_x,v_y,v_z,w_x,w_y,w_z,dim );
            stress_update (sxyptr,c16,c26,c36,c46,c56,c66,z,x,y,dt,u_x,u_y,u_z,v_x,v_y,v_z,w_x,w_y,w_z,dim );
        }
    }
}
#else
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
                                   const dim_t  dim)
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
                const float c11 = cell_coeff_BL      (cc11, z, x, y, dim);
                const float c12 = cell_coeff_BL      (cc12, z, x, y, dim);
                const float c13 = cell_coeff_BL      (cc13, z, x, y, dim);
                const float c14 = cell_coeff_ARTM_BL (cc14, z, x, y, dim);
                const float c15 = cell_coeff_ARTM_BL (cc15, z, x, y, dim);
                const float c16 = cell_coeff_ARTM_BL (cc16, z, x, y, dim);
                const float c22 = cell_coeff_BL      (cc22, z, x, y, dim);
                const float c23 = cell_coeff_BL      (cc23, z, x, y, dim);
                const float c24 = cell_coeff_ARTM_BL (cc24, z, x, y, dim);
                const float c25 = cell_coeff_ARTM_BL (cc25, z, x, y, dim);
                const float c26 = cell_coeff_ARTM_BL (cc26, z, x, y, dim);
                const float c33 = cell_coeff_BL      (cc33, z, x, y, dim);
                const float c34 = cell_coeff_ARTM_BL (cc34, z, x, y, dim);
                const float c35 = cell_coeff_ARTM_BL (cc35, z, x, y, dim);
                const float c36 = cell_coeff_ARTM_BL (cc36, z, x, y, dim);
                const float c44 = cell_coeff_BL      (cc44, z, x, y, dim);
                const float c45 = cell_coeff_ARTM_BL (cc45, z, x, y, dim);
                const float c46 = cell_coeff_ARTM_BL (cc46, z, x, y, dim);
                const float c55 = cell_coeff_BL      (cc55, z, x, y, dim);
                const float c56 = cell_coeff_ARTM_BL (cc56, z, x, y, dim);
                const float c66 = cell_coeff_BL      (cc66, z, x, y, dim);
                
                const float u_x = stencil_X (SX, vxu, dxi, z, x, y, dim);
                const float v_x = stencil_X (SX, vxv, dxi, z, x, y, dim);
                const float w_x = stencil_X (SX, vxw, dxi, z, x, y, dim);
                
                const float u_y = stencil_Y (SY, vyu, dyi, z, x, y, dim);
                const float v_y = stencil_Y (SY, vyv, dyi, z, x, y, dim);
                const float w_y = stencil_Y (SY, vyw, dyi, z, x, y, dim);
                
                const float u_z = stencil_Z (SZ, vzu, dzi, z, x, y, dim);
                const float v_z = stencil_Z (SZ, vzv, dzi, z, x, y, dim);
                const float w_z = stencil_Z (SZ, vzw, dzi, z, x, y, dim);
                
                stress_update (sxxptr,c11,c12,c13,c14,c15,c16,z,x,y,dt,u_x,u_y,u_z,v_x,v_y,v_z,w_x,w_y,w_z,dim );
                stress_update (syyptr,c12,c22,c23,c24,c25,c26,z,x,y,dt,u_x,u_y,u_z,v_x,v_y,v_z,w_x,w_y,w_z,dim );
                stress_update (szzptr,c13,c23,c33,c34,c35,c36,z,x,y,dt,u_x,u_y,u_z,v_x,v_y,v_z,w_x,w_y,w_z,dim );
                stress_update (syzptr,c14,c24,c34,c44,c45,c46,z,x,y,dt,u_x,u_y,u_z,v_x,v_y,v_z,w_x,w_y,w_z,dim );
                stress_update (sxzptr,c15,c25,c35,c45,c55,c56,z,x,y,dt,u_x,u_y,u_z,v_x,v_y,v_z,w_x,w_y,w_z,dim );
                stress_update (sxyptr,c16,c26,c36,c46,c56,c66,z,x,y,dt,u_x,u_y,u_z,v_x,v_y,v_z,w_x,w_y,w_z,dim );
            }
        }
    }
}
#endif

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
                                 const dim_t  dim,
                                 void*        stream)
{
    const int block_dim_x = 32;
    const int block_dim_y = 4;

    dim3 grid_dim( ((nzf-nz0) + block_dim_x-1)/block_dim_x,
                   ((nxf-nx0) + block_dim_y-1)/block_dim_y,
                      1 );
    dim3 block_dim(block_dim_x, block_dim_y, 1);

    cudaStream_t s = (cudaStream_t) stream;

#ifdef SCELL_BL
    const int HALO = 4;
    compute_component_scell_BL_cuda_k<HALO,block_dim_x,block_dim_y><<<grid_dim, block_dim, 0, s>>>
        (sxxptr, syyptr, szzptr, syzptr, sxzptr, sxyptr,
         vxu, vxv, vxw, vyu, vyv, vyw, vzu, vzv, vzw,
         cc11, cc12, cc13, cc14, cc15, cc16,
         cc22, cc23, cc24, cc25, cc26,
         cc33, cc34, cc35, cc36,
         cc44, cc45, cc46,
         cc55, cc56,
         cc66,
         dt, dzi, dxi, dyi, nz0, nzf, nx0, nxf, ny0, nyf, SZ, SX, SY, dim);
#else
    compute_component_scell_BL_cuda_k<<<grid_dim, block_dim, 0, s>>>
        (sxxptr, syyptr, szzptr, syzptr, sxzptr, sxyptr,
         vxu, vxv, vxw, vyu, vyv, vyw, vzu, vzv, vzw,
         cc11, cc12, cc13, cc14, cc15, cc16,
         cc22, cc23, cc24, cc25, cc26,
         cc33, cc34, cc35, cc36,
         cc44, cc45, cc46,
         cc55, cc56,
         cc66,
         dt, dzi, dxi, dyi, nz0, nzf, nx0, nxf, ny0, nyf, SZ, SX, SY, dim);
#endif 
    CUDA_CHECK(cudaGetLastError());
};



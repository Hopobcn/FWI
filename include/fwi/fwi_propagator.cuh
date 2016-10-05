#ifndef _FWI_PROPAGATOR_CUDA_H_
#define _FWI_PROPAGATOR_CUDA_H_

#define CUDA_CHECK(call) { gpu_assert((call), __FILE__, __LINE__); }
inline void gpu_assert(cudaError_t code, const char* file, int line)
{
    if (code != cudaSuccess)
    {
        fprintf(stderr, "CUDA ERROR: %s:%d, ", file, line);
        fprintf(stderr, "code: %d, reason: %s\n", code, cudaGetErrorString(code));
        exit(code);
    }
}

extern "C"
void compute_component_vcell_TL_cuda ( /* COMPLETE ME */ );

extern "C"
void compute_component_vcell_TR_cuda ( /* COMPLETE ME */ );

extern "C"
void compute_component_vcell_BR_cuda ( /* COMPLETE ME */ );

extern "C"
void compute_component_vcell_BL_cuda ( /* COMPLETE ME */ );

extern "C"
void compute_component_scell_TR_cuda ( /* COMPLETE ME */ );

extern "C"
void compute_component_scell_TL_cuda ( /* COMPLETE ME */ );

extern "C"
void compute_component_scell_BR_cuda ( /* COMPLETE ME */ );

extern "C"
void compute_component_scell_BL_cuda ( /* COMPLETE ME */ );

#endif /* end of _FWI_PROPAGATOR_CUDA_H_ definition */

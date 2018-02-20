`FWI` is a great candidate to take advantage of shuffle intrinsics and shared memory.
OpenACC provides the directive `cache` to exploit shared memory, but it lacks a way to exploit CUDA intra-warp intrinsics.

In this step we provide a set of highly optimized kernels in `src/fwi_propagator.cu` file.
Your task consists on adding the necessary glue code to get:
 - The GPU device pointers managed by the OpenACC runtime
 - The CUDA stream allocated by OpenACC

We also provide the necessary modifications in `CMakeLists.txt` for compiling with `nvcc` and linking with `pgcc`.
Just remember to pass `-DUSE_OPENACC=ON -DUSE_CUDA_KERNELS=ON` to cmake.


####In summary:
* Compile `FWI` with `-DUSE_OPENACC=ON -DUSE_CUDA_KERNELS=ON`
* Add `#pragma acc host_data use_device` directives to forward the *device pointers* allocated by OpenACC to our CUDA kernels.
* Pass the current stream to the CUDA kernel (with `acc_get_cuda_stream`).

For instance, for `compute_component_vcell_TL_cuda`:
```c
{
#if !defined(USE_CUDA)
    <... previous OpenACC impl. ...>
#else
    void* stream = acc_get_cuda_stream(phase)

    #pragma acc host_data use_device(szptr, sxptr, syptr, rho, vptr)
    {
        compute_component_vcell_TL_cuda(..., stream);
    }
#endif
};
```

#### Benchmarking

```bash
$ cmake -DCMAKE_C_COMPILER=pgcc -DUSE_OPENMP=OFF -DUSE_OPENACC=ON -DUSE_CUDA_KERNELS=ON ..
$ make irun
[ 27%] Built target fwi-core-cuda
[ 72%] Built target fwi-core
[ 90%] Built target fwi
[100%] outputs will be in /home/ubuntu/FWI/scripts/output/
PROJECT_SOURCE_DIR: /home/ubuntu/FWI
PROJECT_BINARY_DIR: /home/ubuntu/FWI/build/bin
COMPILER_ID:        PGI
---
/home/ubuntu/FWI/build/bin/fwi fwi_profile.txt
---
MPI rank 0 with GPU 0 (1)
Number of frequencies 1
Number of shots 1
Number of gradient iterations 1
Number of test iterations 1
Output directory path: results
FWI Program finished in 7.556227 seconds
[100%] Built target irun
```
We got a 1.07 speedup.
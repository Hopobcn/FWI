In this step, you have to add compute regions around expensive loops in the application. 
We assume you already have some OpenACC knowledge.

In this lab we have used the **kernels** directive extensively but the **loop** directive could also been used.

To reduce time, we will suppose we have profiled the mini-app and we got this results (see Appendix 0):

```bash
Number of frequencies 1
Number of shots 1
Number of gradient iterations 1
Number of test iterations 1
Output directory path: results
FWI Program finished in 1059.575686 seconds

======== CPU profiling result (bottom up):
Time(%)      Time  Name
 19.99%  211.766s  IDX
 16.53%  175.115s  compute_component_scell_TR
 16.28%  172.505s  compute_component_scell_BL
 15.63%  165.634s  compute_component_scell_BR
  9.47%  100.343s  compute_component_scell_TL
  5.87%  62.2217s  compute_component_vcell_TR
  5.83%  61.7416s  compute_component_vcell_BR
  5.17%  54.7615s  compute_component_vcell_BL
  5.15%  54.5714s  compute_component_vcell_TL

======== Data collected at 100Hz frequency
======== Percentage threshold: 1%
```
We can see that `scell` and `vcell` functions dominate the execution time.
`IDX` is the function that linearizes the (i,j,k) triplet into the linear index.
Usually the compiler is smart enough to inline it, but in this execution it didn't.
Since we know that `IDX` is only called inside `scell` and `vcell` functions, we can safely split the `IDX` execution time among `scell` and `vcell` functions.
Therefore we can safely say that `scell` and `vcell` accounts for the 99% of the execution time of the application.

If we take a look at those functions in **src/fwi_propagator.c** we will arrive to this conclusions:
1. They are embarrassingly parallel
2. All `TR`/`TR`/`BL`/`BR` are very similar
3. We have to apply the same parallelization strategy for all *scell* and *vcell* functions.



In this first implementation we are going to use CUDA Unified Memory.
For that purpose we already modified **CMakeLists.txt** to add the `managed` to the `-ta=tesla` openacc target:
```cmake
set(OpenACC_C_FLAGS "${OpenACC_C_FLAGS} -ta=tesla,cuda8.0,cc20,cc35,cc50,cc60,lineinfo,managed")
```



> To facilitate your work we already implemented the majority of `openacc` pragmas leaving `vcell_TL` and `scell_TR` for you to implement.

##### In sumary:
You have to add a `#pragma acc kernels` in the outer-most loop of `src/fwi_propagator.c:166` (`compute_component_vcell_TL`) and `src/fwi_propagator.c:634` (`compute_component_scell_TR`). 
Example for `vcell_TL`:
```c
#pragma acc kernels
for(integer y=ny0; y < nyf; y++)
{
    for(integer x=nx0; x < nxf; x++)
    {
        for(integer z=nz0; z < nzf; z++)
        {
            const real lrho = rho_TL(rho, z, x, y, dimmz, dimmx);

            const real stx  = stencil_X( _SX, sxptr, dxi, z, x, y, dimmz, dimmx);
            const real sty  = stencil_Y( _SY, syptr, dyi, z, x, y, dimmz, dimmx);
            const real stz  = stencil_Z( _SZ, szptr, dzi, z, x, y, dimmz, dimmx);

            vptr[IDX(z,x,y,dimmz,dimmx)] += (stx  + sty  + stz) * dt * lrho;
        }
    }
}
```

Then compile `fwi` (*make sure  to use* **pgcc** *compiler* and to *enable OpenACC* with **-DUSE_OPENACC=ON** !), and observe the compiler output:

```bash
$ mkdir build
$ cd build
$ cmake -DCMAKE_C_COMPILER=pgcc -DUSE_OPENMP=OFF -DUSE_OPENACC=ON ..
$ make
...
[ 11%] Building C object src/CMakeFiles/fwi-core.dir/fwi_propagator.c.o
...
compute_component_vcell_TL:
    166, Generating implicit copy(vptr[:])
    170, Loop carried dependence due to exposed use of vptr[:*] prevents parallelization
         Accelerator scalar kernel generated
         Accelerator kernel generated
         Generating Tesla code
        170, #pragma acc loop seq
        174, #pragma acc loop seq
        180, #pragma acc loop seq
    174, Loop carried dependence due to exposed use of vptr[:*] prevents parallelization
         180, Complex loop carried dependence of vptr-> prevents parallelization
              Loop carried dependence due to exposed use of vptr[:*] prevents parallelization
...
compute_component_scell_TR:
    636, Loop carried dependence due to exposed use of sxzptr[:*],szzptr[:*],syzptr[:*],syyptr[:*],sxyptr[:*],sxxptr[:*] prevents parallelization
         Accelerator scalar kernel generated
         Accelerator kernel generated
         Generating Tesla code
        636, #pragma acc loop seq
        640, #pragma acc loop seq
        646, #pragma acc loop seq
    640, Loop carried dependence due to exposed use of sxzptr[:*],syzptr[:*],szzptr[:*],syyptr[:*],sxyptr[:*],sxxptr[:*] prevents parallelization
         646, Complex loop carried dependence of sxxptr->,syyptr-> prevents parallelization
              Loop carried dependence due to exposed use of sxzptr[:*] prevents parallelization
              Complex loop carried dependence of szzptr-> prevents parallelization
              Loop carried dependence due to exposed use of szzptr[:*],syzptr[:*] prevents parallelization
              Complex loop carried dependence of syzptr-> prevents parallelization
              Loop carried dependence due to exposed use of syyptr[:*] prevents parallelization
              Complex loop carried dependence of sxzptr-> prevents parallelization
              Loop carried dependence due to exposed use of sxyptr[:*] prevents parallelization
              Complex loop carried dependence of sxyptr-> prevents parallelization
              Loop carried dependence due to exposed use of sxxptr[:*] prevents parallelization
...
```

Oops! the compiler detects a dependence and prevents a parallelization (it generates a scalar kernel!). 
Since we know that *vcell_TL* is embarrasingly parallel and there isn't a dependence we have to force the compiler to ignore those dependences and parallelize it.
For that we have to add `#pragma acc loop independent` before each iteration level:
```c
#pragma acc kernels
#pragma acc loop independent
for (integer y=ny0; y < nyf; y++) {
    #pragma acc loop independent
    for (integer x=nx0; x < nxf; x++) {
        #pragma acc loop independent
        for (integer z=nz0; z < nzf; z++) {
```

Then we can compile the application again:
```bash
$ make
Scanning dependencies of target fwi-core
[ 11%] Building C object src/CMakeFiles/fwi-core.dir/fwi_propagator.c.o
...
compute_component_vcell_TL:
    166, Generating implicit copy(vptr[:])
    171, Loop is parallelizable
         Generating Multicore code
        171, #pragma acc loop gang
    176, Loop is parallelizable
         183, Loop is parallelizable
              Accelerator kernel generated
              Generating Tesla code
             171, #pragma acc loop gang /* blockIdx.y */
             176, #pragma acc loop gang, vector(4) /* blockIdx.z threadIdx.y */
             183, #pragma acc loop gang, vector(32) /* blockIdx.x threadIdx.x */
...
compute_component_scell_TR:
    640, Loop is parallelizable
         Generating Multicore code
        640, #pragma acc loop gang
    645, Loop is parallelizable
         652, Loop is parallelizable
              Accelerator kernel generated
              Generating Tesla code
             640, #pragma acc loop gang /* blockIdx.y */
             645, #pragma acc loop gang, vector(4) /* blockIdx.z threadIdx.y */
             652, #pragma acc loop gang, vector(32) /* blockIdx.x threadIdx.x */
...
[  9%] Linking C static library ../lib/libfwi-core.a
[ 33%] Built target fwi-core
[ 38%] Linking C executable ../../bin/fwi-data-generator
[ 42%] Built target fwi-data-generator
[ 47%] Linking C executable ../../bin/fwi-sched-generator
[ 52%] Built target fwi-sched-generator
[ 57%] Linking C executable ../bin/fwi
[ 61%] Built target fwi
[ 76%] Built target Unity
[ 80%] Linking C executable ../bin/fwi-tests
[100%] Built target fwi-tests
```

> Obs: Functions called from a parallel region (`== CUDA kernel`) should be declared a device function. For that OpenACC provides the `#pragma acc routine <type>` directive.
To speed-up the hands-on lab we already provide the code with those pragmas.
We encourage you to check `include/fwi/fwi_propagator.h` to see how it's done.


#### Benchmark

After implementing all *scell* and *vcell* functions we can proceed to measure the execution time:
```bash
$ make irun
[ 62%] Built target fwi-core
[ 87%] Built target fwi
Scanning dependencies of target irun
[100%] outputs will be in /home/ubuntu/FWI/scripts/output/
PROJECT_SOURCE_DIR: /home/ubuntu/FWI
PROJECT_BINARY_DIR: /home/ubuntu/FWI/build/bin
COMPILER_ID:        PGI
---
/home/ubuntu/FWI/build/bin/fwi fwi_schedule.txt
---
MPI rank 0 with GPU 0 (1)
Number of frequencies 1
Number of shots 1
Number of gradient iterations 1
Number of test iterations 1
Output directory path: results
FWI Program finished in 8.664497 seconds
[100%] Built target irun
```
That is, 13.9x faster than the OpenMP execution.

Remember you can see differences with the soluction with `git diff gtc2018-step1-sol`

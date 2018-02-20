OpenACC follows the same semantics as CUDA regarding streams.
By default, kernels and memory copies are executed in the default stream which imposes serialization between kernels and memory transfers.

Since all *vcell** kernels are independent between each other. 
And all *scell** kernels are independent between each other.
We could launch independent kernels to different streams and we could achieve back-to-back execution.
In the worst case, we can begin execution the moment after the previous kernel left some SM empty.
In the best case (when a kernel doesn't occupy all de SMs), multiple kernels could run in concurrently.
Another benefit is that we can perform H2D,D2H & kernel execution in parallel.

OpenACC uses the clause `async` that requires a nonnegative scalar integer as a parameter.
In `include/fwi/fwi_propagator.h::78` we have already prepared an enum that we are going to use as identifiers:
```c
typedef enum {TR, TL, BR, BL} phase_t;
```

We already prepared all `compute_component_scell_*` and `compute_component_vcell*` functions with and additional parameter `phase` that we will use.


##### In sumary:

* Add the `async(phase)` clause to all `acc kernels` in *scell* and *vcell* functions.

For instance for `vcell_TL` (`src/fwi_propagator.c:168`):
```c
#pragma acc kernels ... async(phase)
```
* In `velocity_propagator` and `stress_propagator` functions.
  Modify the last parameter of `compute_component_scell_*` and `compute_component_vcell*` functions to pass `TR`, `TL`, `BR` o `BL` values.


#### Benchmarking

```
# make irun
[ 62%] Built target fwi-core
[ 87%] Built target fwi
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
FWI Program finished in 8.100650 seconds
[100%] Built target irun
```
We got another 1.03x speedup.
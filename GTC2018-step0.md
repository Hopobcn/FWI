The initial step of every paralelization effort is to check whether it makes sense to apply optimizations or not.
In the case of accelerator programming a certain number of questions have to be answered before the first line of code is written.
This includes:
* Understanding the program structure and how data is passed through the call tree
* Profiling the CPU-only version of the application and identifying computationally-intense "hot spots"
* Identify which loop nests dominate the runtime
* Are the loop nests suitable for an accelerator?
* Insuring that the algorithms you are considering for acceleration are safely parallel


Firt. Compile & execute the sequential version of FWI:
```
$ mkdir build
$ cd build
$ cmake -DCMAKE_C_COMPILER=pgcc ..
$ make
```

Then you execute FWI either using `make irun` or `bin/fwi fwi_schedule.txt`. You should see something like this:
```
Number of frequencies 1
Number of shots 1
Number of gradient iterations 1
Number of test iterations 1
Output directory path: results
FWI Program finished in 1059.629833 seconds
```

Then we are going to profile FWI to search for hot spots using `make iprofile` which will call `nvprof --cpu-profiling on` (or call `nvprof` directly from console: `nvprof --cpu-profiling on bin/fwi fwi_schedule.txt`).

```bash
$ cmake -DCMAKE_C_COMPILER=pgcc ..
$ make
$ nvprof --cpu-profiling on --cpu-profiling-percentage-threshold 1 bin/fwi fwi_schedule.txt
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

Now we will recompile enabling OpenMP and execute the application to measure the performance of OpenMP vs the serial implementation:
```bash
Number of frequencies 1
Number of shots 1
Number of gradient iterations 1
Number of test iterations 1
Output directory path: results
FWI Program finished in 120.587904 seconds
```

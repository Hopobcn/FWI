Now discard all non-commited changes and checkout branch `gtc2018-step2`:

Step 2 will profile the application to find possible weaknesses and optimization opportunities. 
We could use *NVIDIA Visual Profiler* for a graphical assestment or `pgprof`/`nvprof` for a command-line visualization. 
For simplicity in this lab we are going to use `nvprof`:
```bash
$ nvprof --dependency-analysis bin/fwi fwi_schedule.txt
==1001== NVPROF is profiling process 1001, command: bin/fwi fwi_schedule.txt
MPI rank 0 with GPU 0 (1)
Number of frequencies 1
Number of shots 1
Number of gradient iterations 1
Number of test iterations 1
Output directory path: results
FWI Program finished in 11.814284 seconds
==1001== Profiling application: bin/fwi fwi_schedule.txt
==1001== Profiling result:
Time(%)      Time     Calls       Avg       Min       Max  Name
 21.65%  2.29883s       450  5.1085ms  456.05us  14.419ms  compute_component_scell_BR_905_gpu
 21.51%  2.28399s       450  5.0755ms  464.31us  14.305ms  compute_component_scell_TR_652_gpu
 20.80%  2.20928s       450  4.9095ms  445.40us  13.852ms  compute_component_scell_BL_1032_gpu
 12.65%  1.34301s       450  2.9845ms  278.01us  8.4036ms  compute_component_scell_TL_778_gpu
  6.10%  647.66ms      1350  479.75us  48.063us  1.3598ms  compute_component_vcell_BR_291_gpu
  5.76%  611.89ms      1350  453.25us  44.063us  1.2821ms  compute_component_vcell_TR_237_gpu
  5.73%  609.01ms      1350  451.12us  44.383us  1.2803ms  compute_component_vcell_BL_345_gpu
  5.66%  601.58ms      1350  445.61us  43.871us  1.2580ms  compute_component_vcell_TL_183_gpu
  0.14%  14.908ms       116  128.52us  128.16us  129.66us  set_array_to_constant_52_gpu

==1001== Unified Memory profiling result:
Device "Tesla K80 (0)"
   Count  Avg Size  Min Size  Max Size  Total Size  Total Time  Name
       1  4.0000KB  4.0000KB  4.0000KB  4.000000KB  3.648000us  Host To Device
       1  4.0000KB  4.0000KB  4.0000KB  4.000000KB  3.552000us  Device To Host
Total CPU Page faults: 1

...

==1001== Dependency Analysis:
==1001== Analysis progress: 100%
Critical path(%)  Critical path  Waiting time  Name
          19.27%      2.298834s           0ns  compute_component_scell_BR_905_gpu
          19.15%      2.283991s           0ns  compute_component_scell_TR_652_gpu
          18.52%      2.209276s           0ns  compute_component_scell_BL_1032_gpu
          11.26%      1.343006s           0ns  compute_component_scell_TL_778_gpu
           5.43%   647.659179ms           0ns  compute_component_vcell_BR_291_gpu
           5.13%   611.886147ms           0ns  compute_component_vcell_TR_237_gpu
           5.10%   608.965553ms           0ns  compute_component_vcell_BL_345_gpu
           5.04%   601.579973ms           0ns  compute_component_vcell_TL_183_gpu
           3.57%   426.061228ms           0ns  cuMemAllocManaged
           2.94%   350.987427ms           0ns  <Other>
           2.14%   254.883255ms           0ns  cuDevicePrimaryCtxRelease
           1.85%   220.329541ms           0ns  cuDevicePrimaryCtxRetain
           0.20%    24.006481ms           0ns  cuMemFree_v2
           0.18%    21.969030ms    10.619674s  cuStreamSynchronize
           ...        ...               ...       ...
...
```
The Critical path is the set of functions which determine the maximum execution time of the application.
Therefore optimization of the critical path should be our first priority.

We can see that *scell* kernels take a good chunck of the critical path.
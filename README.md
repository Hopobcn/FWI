# FWI mini-app

## Parallelization of a Reverse Time Migration (RTM) program using OpenACC/CUDA

Reverse time migration (RTM) modeling is a critical component in the seismic
processing workflow of oil and gas exploration as well as for the understanding
of energy release in subduction-zone earthquakes. With the help of high
precision seismic sensors deployed on the field, it is possible to use the
information gathered during seismic aftershocks and reverse-time migrate them.
This can give scientists a large amount of highly accurate information of the
seismic conditions of the region of interest.

Such analysis is critical after a large earthquake because it can help
scientists know the state of a seismic fault and the probability of subsequent
large aftershocks. As the number of aftershocks sky-rockets after a large
earthquake, the amount of data to analyse grows really fast.  Thus, it is
mandatory to speed up the processing of all that information.

### Build Instructions:

This application uses CMake to discover all dependences and build the application. Therefore the build process follows the typical build process of every cmake application.

```bash
cmake -DCMAKE_C_COMPILER=<foo-compiler> [ -D<OPTION_1>=<yes|no> -D<OPTION_2>=<YES|NO> ... ]  <path-to-project-base-dir>
```

__WARNING:__ *Always* make *out-of-source* builds (don't execute cmake from the project root directory):
```bash
cd ~/<repo-name>/
mkdir build && cd build
cmake <options> ..
make
```


The FWI/CMakeLists.txt has been modified to accept those options:

| CMake Options    | Default Value | Description                           | Observations                             |
| -----------------|:-------------:| ------------------------------------- |------------------------------------------|
| ENABLE_TESTS     | ON            | Build tests                           |                                          |
| PERFORM_IO       | OFF           | Load/Store dataset from disc          | Should be OFF when measuring performance |
| IO_STATS         | OFF           | Log fwrite/fread performance          |                                          |
| USE_MPI          | OFF           | Enable MPI compilation                |                                          |
| USE_OPENMP       | OFF           | Enable OpenMP compilation             | Either OpenMP or OpenACC must be enabled  not both |
| USE_OPENACC      | OFF           | Enable OpenACC compilation            | OpenACC requires the PGI 16.5 compiler   |
| USE_CUDA_KERNELS | OFF           | Enable CUDA kernels back-end          | Requires OpenACC to be enabled           |
| PROFILE          | OFF           | Add profile information to the binary |                                          |

We provide three scripts `environment_icc.sh`, `environment_gcc.sh` and `environment_pgi.sh` that will load the required modules to compile the application with the supported compilers.
But if you feel lucky you can load any other modules available in Minotauro. `module avail` will display all applications/compilers/tools available in this cluster which can be loaded with `module load <name>` and unloaded with `module unload <name>`.
For more information read the [Minotauro User Gudie](http://www.bsc.es/user-support/mt.php)

We also provide a matrix of the three compilers supported (`icc >= 16.0.2`, `gcc >= 4.9.3` and `pgcc >= 16.5`) with all CMake valid options:

| CMake Options    | ICC    | GCC    | PGI 16.5 |
| -----------------|:------:|:------:|:------:|
| ENABLE_TESTS     | YES    | YES    | YES    |
| PERFORM_IO       | YES    | YES    | YES    |
| IO_STATS         | YES    | YES    | YES    |
| USE_MPI          | YES    | YES    | YES    |
| USE_OPENMP       | YES    | YES    | YES*   |
| USE_OPENACC      | NO     | NO     | YES*   |
| USE_CUDA_KERNELS | NO     | NO     | YES    |
| PROFILE          | YES    | YES    | YES    |

*The code is prepared to use OpenMP parallelization or OpenACC acceleration, not both at the same time, so please use only one option at build time.

CMake also provides a set of BUILD_TYPES. In our case we will use `Release` for performance tests and `Debug` for debugging the application (the default is set to `Release`). Example:
```bash
cmake -DCMAKE_C_COMPILER=<foo-compiler> -DCMAKE_BUILD_TYPE=Debug [ -D<OPTION_1>=<yes|no> -D<OPTION_2>=<YES|NO> ]  <path-to-project-base-dir>
make
```

#### Examples:

Building FWI sequential with ICC:

```bash
source scripts/environment_icc.sh
cd build
cmake -DCMAKE_C_COMPILER=icc ..
make
```

Building FWI sequential with GCC with profiling support:

```bash
source scripts/environment_gcc.sh
cd build
cmake -DCMAKE_C_COMPILER=gcc -DCMAKE_BUILD_TYPE=Release -DPROFILE=ON ..
make
```

Building FWI with OpenMP support with PGI and execute tests:
```bash
source scripts/environment_pgi.sh
cd build
cmake -DCMAKE_C_COMPILER=pgcc -DCMAKE_BUILD_TYPE=Release -DUSE_OPENMP=YES ..
make utest
```

Building FWI with OpenACC+CUDA kernels with PGI and execute tests:
```bash
source scripts/environment_pgi.sh
cd build
cmake -DCMAKE_C_COMPILER=pgcc -DCMAKE_BUILD_TYPE=Release -DUSE_OPENMP=NO -DUSE_OPENACC=YES -DUSE_CUDA_KERNELS=YES ..
make utest
```

#### Running Instructions:

To facilitate your work launching jobs, we added a set of 'targets' that launch some SLURM scripts in the queue system.

|  Makefile target  | acction                                          | example       |
| -----------------|:------------------------------------------------|:--------------|
| run-seq          | Launches script/jobscript_run.sequential.slurm   | `make run-seq` |
| run-openmp       | Launches script/jobscript_run.openmp.slurm       | `make run-openmp` |
| run-openacc      | Launches script/jobscript_run.openacc.slurm      | `make run-openacc` |

All executions **should** be performed in compute nodes. You can create/modify any script in `scripts` folder but **do not** modify the execution wall time.

   
#### Profiling Instructions:

To help you, CMakeLists.txt has been modified to include `-pg` (gcc), `-p` (Intel) or `-Mprof` (PGI) when `-DPROFILE=ON` is provided:
```bash
source scripts/environment_gcc.sh
cd build
cmake -DCMAKE_C_COMPILER=gcc -DCMAKE_BUILD_TYPE=Release -DPROFILE=ON ..
make
<--- run application (in compute node not login node!) --->
<--- run gprof                                         --->
```

### References

[Minotauro User Gudie](http://www.bsc.es/user-support/mt.php)

[The OpenACC Application Programmin Interface V2.5](http://www.openacc.org/sites/default/files/OpenACC_2pt5.pdf)

[CUDA C Programming Guide](http://docs.nvidia.com/cuda/cuda-c-programming-guide)

[CMAKE Documentation](https://cmake.org/cmake/help/v3.6/)

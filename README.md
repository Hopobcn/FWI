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

## [GTC2017: Best GPU Code Practices Combining OpenACC, CUDA and OmpSs - Lab Instructions](GTC2017.md)

### How to get the code:

Pull the repo using `git clone`:
```bash
git clone https://github.com/Hopobcn/FWI.git
```
Download all pre-requisites using `git submodules`:
```bash
git submodule update --init --recursive
```

### Build Instructions:

This application uses CMake to discover all dependences and build the application. Then the build process follows the typical build process of every cmake application:

```bash
cmake -DCMAKE_C_COMPILER=<foo-compiler> -DCMAKE_BUILD_TYPE=<Release|Debug> [ -D<OPTION_1>=<yes|no> -D<OPTION_2>=<YES|NO> ... ]  <path-to-project-base-dir>
```

__WARNING:__ *Always* make *out-of-source* builds (don't execute cmake from the project root directory):
```bash
cd ~/FWI/
mkdir build && cd build
cmake <options> ..
make
```

#### Build Options:

The FWI/CMakeLists.txt has been modified to accept those options:

| CMake Options    | Default Value | Description                           | Observations                             |
| -----------------|:-------------:| ------------------------------------- |------------------------------------------|
| ENABLE_TESTS     | ON            | Build tests                           |                                          |
| PERFORM_IO       | OFF           | Load/Store dataset from disc          | Should be OFF when measuring performance |
| IO_STATS         | OFF           | Log fwrite/fread performance          |                                          |
| USE_MPI          | OFF           | Enable MPI compilation                |                                          |
| USE_OPENMP       | OFF           | Enable OpenMP compilation             | Either OpenMP or OpenACC must be enabled  not both |
| USE_OPENACC      | OFF           | Enable OpenACC compilation            | OpenACC requires the PGI +16.5 compiler  |
| USE_CUDA_KERNELS | OFF           | Enable CUDA kernels back-end          | Requires OpenACC to be enabled           |
| PROFILE          | OFF           | Add profile information to the binary |                                          |


#### How to execute FWI:

Usage:
```bash
bin/fwi <params-file> <frequency-file>
```
Example:
```bash
bin/fwi ../data/fwi_params.txt ../data/fwi_frequencies.profile.txt
```

#### CPU Profiling Instructions:

To profile the CPU execution, use `-DPROFILE=ON` to include `-pg` (gcc), `-p` (Intel) or `-Mprof` (PGI) automatically:
```bash
cmake -DCMAKE_C_COMPILER=gcc -DCMAKE_BUILD_TYPE=Release -DPROFILE=ON ..
make
<--- run application --->
<--- run gprof       --->
```

#### Examples:

Build FWI sequential with ICC:

```bash
source scripts/environment_icc.sh
cd build
cmake -DCMAKE_C_COMPILER=icc ..
make
```

Build FWI sequential with GCC with profiling support:

```bash
source scripts/environment_gcc.sh
cd build
cmake -DCMAKE_C_COMPILER=gcc -DCMAKE_BUILD_TYPE=Release -DPROFILE=ON ..
make
```

Build FWI with OpenMP support with PGI and execute tests:

```bash
source scripts/environment_pgi.sh
cd build
cmake -DCMAKE_C_COMPILER=pgcc -DCMAKE_BUILD_TYPE=Release -DUSE_OPENMP=YES ..
make utest
```

Build FWI with OpenACC with PGI and execute tests:

```bash
source scripts/environment_pgi.sh
cd build
cmake -DCMAKE_C_COMPILER=pgcc -DCMAKE_BUILD_TYPE=Release -DUSE_OPENMP=NO -DUSE_OPENACC=YES ..
make utest
```

Build FWI with OpenACC+CUDA kernels with PGI and execute tests:

```bash
source scripts/environment_pgi.sh
cd build
cmake -DCMAKE_C_COMPILER=pgcc -DCMAKE_BUILD_TYPE=Release -DUSE_OPENMP=NO -DUSE_OPENACC=YES -DUSE_CUDA_KERNELS=YES ..
make utest
```

### Minotauro:

#### Build FWI in Minotauro:

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


#### Run FWI in Minotauro:

To facilitate your work launching jobs, we added a set of 'targets' that launch some SLURM scripts in the queue system.

|  Makefile target  | acction                                          | example       |
| -----------------|:------------------------------------------------|:--------------|
| run-seq          | Launches script/jobscript_run.sequential.slurm   | `make run-seq` |
| run-openmp       | Launches script/jobscript_run.openmp.slurm       | `make run-openmp` |
| run-openacc      | Launches script/jobscript_run.openacc.slurm      | `make run-openacc` |

All executions **should** be performed in compute nodes. You can create/modify any script in `scripts` folder but **do not** modify the execution wall time.


### References

[Minotauro User Gudie](http://www.bsc.es/user-support/mt.php)

[The OpenACC Application Programmin Interface V2.5](http://www.openacc.org/sites/default/files/OpenACC_2pt5.pdf)

[CUDA C Programming Guide](http://docs.nvidia.com/cuda/cuda-c-programming-guide)

[CMAKE Documentation](https://cmake.org/cmake/help/v3.6/)

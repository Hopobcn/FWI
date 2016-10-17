# BSC/UPC Hackathon

## GPU Part - Parallelization of a Reverse Time Migraton (RTM) program using OpenACC/CUDA

Reverse time migration (RTM) modeling is a critical component in the seismic processing workflow of oil and gas exploration.

In this Hackathon you will parallelize a mini-app called FWI using OpenACC and CUDA.

In order to get the maximum points, ten questions will be presented in this document that you will have to answer and perform the steps asked. Each question/step has a dificulty associated, from the easiest one (question 1) to  the hardest (question 10).

Please follow that order. Each step has to be documented using `git`. Every time you consider you have a good answer for question 'N', then proceed to do a`git commit -m "Answer N"` & `git push` in `master` branch.

### How to get started:
`clone` your repository into your local machine:
```bash
git clone https://hackathon.hca.bsc.es/hackathon/<username>.git
git submodule update --init --recursive
```
Then upload the whole directory to Minotauro using `scp` or `rsync`.
There is no connectivity inside Minotauro to the outside world so you will have to copy the entire directory from Minotauro to your laptop to be able to push to github.

Example using `rsync`:

Local to Minotauro:
```bash
rsync -azP --delete --exclude 'build' /home/foo/hackathon/FWI/ USER@mt2.bsc:~/hackathon/FWI
```
Minotauro to Local:
```bash
rsync -azP --delete --exclude 'build' USER@mt2.bsc:~/hackathon/FWI /home/foo/hackathon
```

### Build Instructions:

This application uses CMake to discover all dependences and build the application. Therefore the build process follows the typical build process of every cmake application.

```bash
cmake -DCMAKE_C_COMPILER=<foo-compiler> [ -D<OPTION_1>=<yes|no> -D<OPTION_2>=<YES|NO> ... ]  <path-to-project-base-dir>
```

__WARNING:__ *Always* make *out-of-source* builds (don't execute cmake from the project root directory):
```bash
cd ~/hackathon/FWI
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

To facilitate your job we provide three scripts `environment_icc.sh`, `environment_gcc.sh` and `environment_pgi.sh` that will load the required modules to compile the application with the supported compilers.
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

Building FWI sequential with ICC 16.0.2:

```bash
source scripts/environment_icc.sh
cd build
cmake -DCMAKE_C_COMPILER=icc ..
make
```

Building FWI sequential with GCC 6.1.0 with profiling support:

```bash
source scripts/environment_gcc.sh
cd build
cmake -DCMAKE_C_COMPILER=gcc -DCMAKE_BUILD_TYPE=Release -DPROFILE=ON ..
make
```

Building FWI with OpenMP support with PGI 16.5 and execute tests:
```bash
source scripts/environment_pgi.sh
cd build
cmake -DCMAKE_C_COMPILER=pgcc -DCMAKE_BUILD_TYPE=Release -DUSE_OPENMP=YES ..
make utest
```

Building FWI with OpenACC+CUDA kernels with PGI 16.5 and execute tests:
```bash
source scripts/environment_pgi.sh
cd build
cmake -DCMAKE_C_COMPILER=pgcc -DCMAKE_BUILD_TYPE=Release -DUSE_OPENMP=NO -DUSE_OPENACC=YES -DUSE_CUDA_KERNELS=YES ..
make utest
```


### Questions/Steps:

Remember to **document** each step using `git`.

All executions **should** be performed in compute nodes. We provide some scripts under the folder `scripts` that can be copied and modified by you to do whatever you want.
**Do not** modify the execution wall time.

1. **Profile the sequential program and give a report of the most time consuming parts of this application. Also tell which of those parts should be ported to the GPU and why.**


    You have the liberty to choose which profiler/tool to use. Some possible options would be:

    * GNU profiler (gprof)
    * NVIDIA profiler (nvprof) (it supports CPU sampling)
    * BSCTOOLS
    * Any other

    CMakeLists.txt has been modified to include `-pg` (gcc), `-p` (Intel) or `-Mprof` (PGI) when `-DPROFILE=ON` is provided:
    ```bash
    source scripts/environment_gcc.sh
    cd build
    cmake -DCMAKE_C_COMPILER=gcc -DCMAKE_BUILD_TYPE=Release -DPROFILE=ON ..
    make
    <--- run application (in compute node not login node!) --->
    <--- run gprof                                         --->
    ```

    We also added support for launching SLURM jobs  (`scripts/job_*.slurm`) with cmake/make (if you want to tweak this go to `main/CMakeLists.txt`):
    ```bash
    # to launch script/jobscript_<COMPILER>_run.sequential.slurm execute:
    make run-seq
    <-- or -->
    # to launch script/jobscript_<COMPILER>_run.openmp.slurm execute:
    make run-openmp
    <-- or -->
    # to launch script/jobscript_<COMPILER>_run.openacc.slurm execute:
    make run-openacc
     ```
    Otherwise use `mnsubmit` to launch your jobs.

2. **Parallelize with OpenACC with 1 GPU**

    From now on, use the PGI 16.5 compiler to work with OpenACC since it's the only fully supported compiler with OpenACC 2.5. (`gcc` >= 6.1.0 has OpenACC 2.0 support but has not been tested)
    We recommend reading the document "[The OpenACC Application Programmin Interface V2.5](www.openacc.org/sites/default/files/OpenACC_2pt5.pdf)" information about the OpenACC Spec 2.5.
    Again, remember to document all the code that you add.

    Add preprocessor guards like:
    ```
    #if defined(_OPENACC)
    #pragma acc <---etc-->
    #elif defined(_OPENMP)
    #pragma omp for <--etc-->
    #endif
    ```
    This way you can recompile your application enabling/disabling this new functionality without breaking the previous implementation (seq/openmp).

3. **Test your OpenACC implementation**

    The serial implementation comes with a set of tests that checks the correct execution of some parts of the program:
    ```bash
    make utest
    ```
    In step (2) you have added some OpenACC pragmas, and now you have to modify the tests in order to make them pass.
    First check that sequential execution passes all tests and then go ahead with your OpenACC implementation.

    In later steps, if you make changes to the code that fall out of the socope of the original tests, you should implement your own tests to be sure that the program runs as expected.

4. **Optimize the OpenACC performance**

    After you know your implementation is correct (remember, Speedup of an incorrect code is ZERO), proceed with optimizing the performance of your kernels.

    For reference purposes this is the execution time of our OpenACC implementation with 1 K80 and freq 10.0 Hz:
    19.7872 s


5. **Implement a Multi-GPU implementation**

    Once you have an optimized single-GPU OpenACC implementation you can proceed to implement a Multi-GPU implementation. Each Minotauro node has 4 K80.
    We recommend using MPI+OpenACC to get a multi-gpu execution.

    Put preprocessor guards as much as possible to prevent braking the implementations that do not use MPI/OpenMP.

6. **Test your Multi-GPU implementation**

    Explain how you checked the correctness.

7. **Optimize the Multi-GPU implementation**

    Optimize your implementation and report execution times of your solution.

8. **Implement velocity (vcell) kernels with CUDA + Tests**

    Instead of writing a new FWI-CUDA application from scratch, we recommend using the OpenACC directive *host_data* that makes the address of a device data available on the host.
    This allows OpenACC/CUDA interoperability, making possible to implement CUDA Kernels while still using OpenACC pragmas in the rest of your code.
    More information in OpenACC 2.5 spec and [openacc-interoperability](github.com/jefflarkin/openacc-interoperability)

    So implement your hand-made kernels and test them for only 1 GPU.

9. **Implement stress (scell) kernels with CUDA + Tests**

    Implement the 'scell' kernels and thest them for only 1 GPU.

10. **Optimize all CUDA kernels**

    Optimize all CUDA kernels that you implemented in case you feel they have to be improved somehow.

    Report kernel times.

    For reference purposes this is the execution time of our OpenACC+CUDA implementation with 1 K80 and freq 10.0 Hz:
    20.4410 s


### References

[Minotauro User Gudie](http://www.bsc.es/user-support/mt.php)

[The OpenACC Application Programmin Interface V2.5](http://www.openacc.org/sites/default/files/OpenACC_2pt5.pdf)

[CUDA C Programming Guide](http://docs.nvidia.com/cuda/cuda-c-programming-guide)

[CMAKE Documentation](https://cmake.org/cmake/help/v3.6/)

# BSC/UPC Hackathon

## GPU Part - Parallelization of a Reverse Time Migraton (RTM) program using OpenACC/CUDA

Reverse time migration (RTM) modeling is a critical component in the seismic processing workflow of oil and gas exploration.

In this Hackathon you will parallelize a mini-app called FWI using OpenACC and CUDA.

In order to get the maximum points, ten questions will be presented in this document that you will have to answer and perform the steps asked. Each question/step has a dificulty associated, from the easiest one (question 1) to  the hardest (question 10).

Please follow that order. Each step has to be documented using `git`. Every time you consider you have a good answer for question 'N', then proceed to do a`git commit -m "Answer N"` & `git push` in `master` branch.

### Build Instructions:

This application uses CMake to discover all dependences and build the application. Therefore the build process follows the typical build process of every cmake application.

```bash
cmake -DCMAKE_C_COMPILER=<foo-compiler> [ -D<OPTION_1>=<yes|no> -D<OPTION_2>=<YES|NO> ]  <path-to-project-base-dir>
```

__WARNING:__ *Always* make *out-of-source* builds (don't execute cmake from the project root directory):
```bash
cd FWI-PROJECT-DIR
mkdir build && cd build
cmake <options> ..
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

CMake also provides a set of BUILD_TYPES. In our case we will use `Release` for performance tests and `Debug` for debugging the application (default is set to `Release`). Example:
```bash
cmake -DCMAKE_C_COMPILER=<foo-compiler> -DCMAKE_BUILD_TYPE=Debug [ -D<OPTION_1>=<yes|no> -D<OPTION_2>=<YES|NO> ]  <path-to-project-base-dir>
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

1. Profile the sequential program and give a report of the most time consuming parts of this application. Also tell which of those parts should be ported to the GPU and why.

    You have the liberty to choose which profiler/tool to use.

2. Parallelize with OpenACC with 1 GPU

    From now on, use the PGI 16.5 compiler to work with OpenACC since it's the only fully supported compiler with OpenACC 2.5. (`gcc` >= 6.1.0 has OpenACC 2.0 support but has not been tested)
    We recommend reading the document "[The OpenACC Application Programmin Interface V2.5](www.openacc.org/sites/default/files/OpenACC_2pt5.pdf)" information about the OpenACC Spec 2.5.
    Again, remember to document all the code that you add.

3. Test your OpenACC implementation

    The serial implementation comes with a set of tests that check the correct execution of some parts of the program:
    ```bash
    make utest
    ```
    In step (2) you have added some OpenACC pragmas, and now you have to modify the tests in order to make them pass.
    First check that sequential execution passes all tests and then go ahead with your OpenACC implementation.

    In later steps, if you make changes to the code that fall out of the socope of the original tests, you should implement your own tests to be sure that the program runs as expected.

4. Optimize the OpenACC performance

    After you know your implementation is correct (remember, Speedup of an incorrect code is ZERO), proceed with optimizing the performance of your kernels.

    For reference purposes this is the execution time of our OpenACC implementation with 1 K80 and freq 10.0 Hz:
    19.7872 s
    With `nvprof`:

    ==31606== Profiling application: ../../build/bin/fwi ../../data/fwi_params.txt ../../data/fwi_frequencies.medium.txt
    ==31606== Profiling result:
    Time(%)      Time     Calls       Avg       Min       Max  Name
     22.36%  6.18451s      1500  4.1230ms  421.53us  12.123ms  compute_component_scell_TR_707_gpu
     21.85%  6.04346s      1500  4.0290ms  414.20us  12.329ms  compute_component_scell_BL_1189_gpu
     21.25%  5.87930s      1500  3.9195ms  379.23us  11.288ms  compute_component_scell_BR_1028_gpu
     12.45%  3.44318s      1500  2.2955ms  225.47us  6.4109ms  compute_component_scell_TL_867_gpu
      5.91%  1.63478s      4500  363.28us  39.935us  1.1306ms  compute_component_vcell_BR_306_gpu
      5.41%  1.49524s      4500  332.27us  36.543us  978.70us  compute_component_vcell_TR_234_gpu
      5.38%  1.48904s      4500  330.90us  36.320us  921.74us  compute_component_vcell_BL_378_gpu
      5.38%  1.48700s      4500  330.44us  36.767us  1.1637ms  compute_component_vcell_TL_162_gpu
      0.01%  3.1546ms        34  92.781us  92.671us  93.407us  set_array_to_random_real_33_gpu
      0.01%  2.2274ms        24  92.809us  92.670us  93.374us  set_array_to_constant_45_gpu
      0.00%  242.01us        59  4.1010us  3.2320us  4.4150us  [CUDA memcpy HtoD]
    (...)


5. Implement a Multi-GPU implementation.

    Once you have an optimized single-GPU OpenACC implementation you can proceed to implement a Multi-GPU implementation. Each Minotauro node has 4 K80.
    We recommend using MPI+OpenACC to get a multi-gpu execution.

6. Test your Multi-GPU implementation
    Explain how you checked the correctness.

7. Optimize the Multi-GPU implementation
    Optimize your implementation and report execution times of your solution.

8. Implement velocity (vcell) Kernels with CUDA + Tests

    Instead of writing a new FWI-CUDA application from scratch, we recommend using the OpenACC directive "host_data" that makes the address of a device data available on the host.
    This allows OpenACC/CUDA interoperability, making possible to implement CUDA Kernels while still using OpenACC pragmas in the rest of your code.
    More information in OpenACC 2.5 spec and [openacc-interoperability](github.com/jefflarkin/openacc-interoperability)

    So implement your hand-made kernels and test them for only 1 GPU.

9. Implement stress (scell) kernels with CUDA + Tests

    Implement the 'scell' kernels and thest them for only 1 GPU.

10. Optimize all CUDA kernels

    Optimize all CUDA kernels that you implemented in case you feel they have to be improved somehow.

    Report kernel times.

     For reference purposes this is the execution time of our OpenACC+CUDA implementation with 1 K80 and freq 10.0 Hz:
    20.4410 s
    With `nvprof`:

    ==31986== Profiling application: ../../build/bin/fwi ../../data/fwi_params.txt ../../data/fwi_frequencies.medium.txt
    ==31986== Profiling result:
    Time(%)      Time     Calls       Avg       Min       Max  Name
     30.97%  11.3465s      1500  7.5644ms  819.51us  15.005ms  void compute_component_scell_TR_cuda_k
     19.72%  7.22522s      1500  4.8168ms  299.07us  9.6079ms  void compute_component_scell_BL_cuda_k
     19.27%  7.05922s      1500  4.7062ms  386.04us  13.415ms  void compute_component_scell_BR_cuda_k
      9.45%  3.46312s      1500  2.3087ms  197.98us  6.3949ms  void compute_component_scell_TL_cuda_k
      5.84%  2.13947s      4500  475.44us  50.176us  1.3116ms  void compute_component_vcell_BR_cuda_k
      5.14%  1.88446s      4500  418.77us  42.271us  1.2819ms  void compute_component_vcell_TL_cuda_k
      5.03%  1.84235s      4500  409.41us  43.743us  1.2252ms  void compute_component_vcell_TR_cuda_k
      4.55%  1.66689s      4500  370.42us  40.032us  1.0153ms  void compute_component_vcell_BL_cuda_k
      0.01%  3.1551ms        34  92.796us  92.671us  93.566us  set_array_to_random_real_33_gpu
      0.01%  2.2275ms        24  92.811us  92.702us  93.214us  set_array_to_constant_45_gpu
      0.00%  242.68us        59  4.1130us  3.2640us  4.5760us  [CUDA memcpy HtoD]


### References

[Minotauro User Gudie](http://www.bsc.es/user-support/mt.php)
[The OpenACC Application Programmin Interface V2.5](www.openacc.org/sites/default/files/OpenACC_2pt5.pdf)
[CUDA C Programming Guide](http://docs.nvidia.com/cuda/cuda-c-programming-guide)
[CMAKE Documentation](https://cmake.org/cmake/help/v3.6/)

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

# [GTC 2018 Instructions](GTC2018.md)

## Getting Started

### Prerequisites:

General prerequisites:
* CMake 3.8 or later
* C compiler (Tested with `gcc` and `pgcc`, `clang` and `icc` should also work with few modifications)

OpenACC prerequisites:
* PGI 17.4 or later

CUDA prerequisites:
* OpenACC prerequisites
* CUDA 8 or later (`nvcc`)

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


__WARNING:__ Advice: *Always* make *out-of-source* builds (don't execute cmake from the project root directory):

Create a `build` directory:
```bash
mkdir build
cd build
```
Execute CMake to generate all the Makefiles required to compile the sequential implementation. And type `make` to compile the application.
```
cmake ..
make
```

#### Build Options:

| CMake Options    | Default Value | Description                           | Observations                             |
| -----------------|:-------------:| ------------------------------------- |------------------------------------------|
| ENABLE_TESTS     | OFF           | Build tests                           | Requires git submodule `Unity`                                         |
| USE_MPI          | OFF           | Enable MPI compilation                |                                          |
| USE_OPENMP       | OFF           | Enable OpenMP compilation             | Either OpenMP or OpenACC must be enabled  not both |
| USE_OPENACC      | OFF           | Enable OpenACC compilation            | Requires compiler with OpenACC 2.5 or above  |
| USE_CUDA_KERNELS | OFF           | Enable CUDA kernels back-end          | Requires OpenACC to be enabled           |
| PROFILE          | OFF           | Add profile information to the binary |                                          |
| PERFORM_IO       | OFF           | Load/Store dataset from disc          | Should be OFF when measuring performance |
| IO_STATS         | OFF           | Log fwrite/fread performance          |                                          |


### Some examples:

> OBS: sometimes, when changing options, CMake may complain or even refuse to generate the Makefiles. One safe solution is to remove all `build/*' contents and try again

1. Sequential
```
cmake ..
make
```
2. OpenMP
```
cmake -DUSE_OPENMP=ON ..
make
```

3. OpenACC
```
cmake -DCMAKE_C_COMPILER=pgcc -DUSE_OPENACC=ON ..
make
```

4. OpenACC + CUDA kernels
```
cmake -DCMAKE_C_COMPILER=pgcc -DUSE_OPENACC=ON -DUSE_CUDA_KERNELS=ON ..
make
```


### How to execute FWI:

The `fwi` binary (located in the `build/bin` directory) depends on a single file: `fwi_schedule.txt` (located in `data` directory).
The application also expects `FWIDIR` env var to be set to the project root directory:

```bash
export FWIDIR=/path/to/FWI/
bin/fwi fwi_schedule.txt
```

The `fwi_schedule.txt` is generated using the `fwi-sched-generator` which depends on `fwi_params.txt` and `fwi_frequencies.txt` files.
Don't modify `fwi_schedule.txt` directly. 
If you wish to modify it, you should do so by modifiying `params` and `frequiencies` files and execute the generator:
```bash
bin/fwi-sched-generator fwi_params.txt fwi_frequencies.txt
```

#### CPU Profiling Instructions:

To profile the CPU execution, use `-DPROFILE=ON` to include `-pg` (gcc), `-p` (Intel) or `-Mprof` (PGI) automatically:
```bash
cmake -DCMAKE_C_COMPILER=gcc -DCMAKE_BUILD_TYPE=Release -DPROFILE=ON ..
make
<--- run application --->
<--- run gprof       --->
```

## Authors
* **Samuel Rodríguez** Original work - [github.com/srodrb/FWI](https://github.com/srodrb)
* **Pau Farré** Acceleration using GPUs/OpenACC+CUDA

See also the list of [contributors](https://github.com/hopobcn/FWI/contributors) who participated in this project.

## License
This project is licensed under the BSD-3 License - see the [LICENSE.md](LICENSE.md) for details.
## References

[The OpenACC Application Programmin Interface V2.5](http://www.openacc.org/sites/default/files/OpenACC_2pt5.pdf)

[CUDA C Programming Guide](http://docs.nvidia.com/cuda/cuda-c-programming-guide)

[CMAKE Documentation](https://cmake.org/cmake/help/v3.8/)

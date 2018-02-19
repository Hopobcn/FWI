Unified Memory (UM) can be very inefficient in older GPU generations and an experienced programmer with detailed knowledge of the application algorithm could outperform the Unified Memory.

As we have seen in Step 2, `FWI` doesn't specially suffer from too much UM traffic.
But as an example, in this step we'll migrate from UM to a movement of data using OpenACC directives.

OpenACC offers two sets of directives to to move data from host to device depending on the `scope`:
* `#pragma acc data` clause can be used to delcare a `scope` where the data resides in the GPU.
* `#pragma enter data` / `#pragma exit data` for variables that don't have a clear scope (for instance, specially useful for C++ Class constructors & desctructors).

In our case, allocations and deallocations happens in different scopes.
Then, we use `#pragma acc enter data create` and `#pragma acc exit data delete` to increase the locality in the GPU.

##### In sumary:

* Remove the `managed` flag from the `OpenACC_C_FLAGS` in `CMakeLists.txt` (already done).
* Annotate, using OpenACC directives, which arrays will be use on the GPU:

In `alloc_memory_shot` function (`src/fwi_kernel.c`), *after* allocations (`malloc`), declare all missing arrays for `coeff_t` (we already filled `v_t` and `s_t` ones).
Example:

```c
coeff_t cc = *c;
#pragma acc enter data create(cc)
#pragma acc enter data create(cc.c11[:cellsInVolume])
#pragma acc enter data create(cc.c12[:cellsInVolume])
#pragma acc enter data create(/* COMPLETEME */)
#pragma acc enter data create(/* COMPLETEME */)
... // continue with c13,c14,c15,c16, c22,c23,c24,c25,c26, c33,c34,c35,c36, c44,c45,c46, c55,c56, c66
```

In `free_memory_shot` function, *before* all dealocations (`free`) we should first, deallocate the GPU memory with:

```c
#pragma acc wait

#pragma acc exit data delete(c->c11)
#pragma acc exit data delete(c->c12)
#pragma acc exit data delete(/* COMPLETEME */)
#pragma acc exit data delete(/* COMPLETEME */)
...
#pragma acc exit data delete(c)
...
```

* On each parallel region (`acc kerenels`) you have inform the compiler that all the arrays are already on the GPU.
 You can do that with the `present` clause:

For instace in *vcell_TL*:
```c
#pragma acc kernels present(szptr, sxptr, syptr, rho, vptr)
```

And for *scell_TR* kernels:
```c
#pragma acc kernels present(sxxptr, syyptr, szzptr, syzptr, sxzptr, sxyptr) \
                    present(vxu, vxv, vxw)  \
                    present(vyu, vyv, vyw)  \
                    present(vzu, vzv, vzw)  \
                    present(cc11, cc12, cc13, cc14, cc15, cc16) \
                    present(cc22, cc23, cc24, cc25, cc26) \
                    present(cc33, cc34, cc35, cc36) \
                    present(cc44, cc45, cc46) \
                    present(cc55, cc56) \
                    present(cc66)
```

#### Benchmarking

Now we can run the application again:
```bash
$ make irun
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
FWI Program finished in 8.401898 seconds
[100%] Built target irun
```
We got a minor speedup of 1.02x. 

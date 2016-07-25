This is the shared memory implementation of the FWI mock-up code.

Compilation of the code.
---------------------------
The compilation process is simple, just modify the Makefile according to your preferences and type "make".

At the beginning of the Makefile there are three compilation flags:
COMPILER: by default, the Intel C compiler is used but it is possible to use the GNU compiler as well.
PERFORM_IO: enables / disables the input output. When disabled, the application does not save velocity fields to disk. This may be useful to measure the performance of the stencil computation.
DEBUG: prints some additional information.

Input parameters configuration.
-------------------------------
The application needs two files to run; a fwi_frequencies.txt file and a fwi_params.txt. These files are located on ../SetupParams directory, for further information on how to define them properly look at ../SetupParams/README.txt

If the program is performing I/O (i.e. it was compiled with PERFORM_IO=YES flag) the program may look for a initial velocity model located at ../InputModels. In order to provide a synthetic model for the input, there is a utility called ModelGenerator on this
directory. Just type "make input" and it will generate all velocity input models needed to run the program.


Execution of the code:
------------------------
Just type "make run".



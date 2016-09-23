# Load with cmd: source <this-file-name>

export MPI_HOME=/opt/pgi/linux86-64/2016/mpi/openmpi-1.10.2
export CUDA_HOME=/usr/local/cuda
export ACC_HOME=/opt/pgi/linux86-64/16.5
export PATH=$ACC_HOME/bin:$MPI_HOME/bin:$CUDA_HOME/bin:$PATH
export LD_LIBRARY_PATH=$ACC_HOME/lib:$MPI_HOME/lib:$CUDA_HOME/lib64:$LD_LIBRARY_PATH
# For OpenACC runtime prfiling with github.com/Hopobcn/openacc-nvtx library
#export ACC_PROFLIB=/home/users/pfarre/scratch-local/CASE/openacc-nvtx/liboaccnvtx.so

# Load with cmd: source <this-file-name>

export MPI_HOME=/opt/pgi/linux86-64/2016/mpi/openmpi-1.10.2
export CUDA_HOME=/usr/local/cuda
export PATH=/opt/pgi/linux86-64/16.5/bin:$MPI_HOME/bin:$CUDA_HOME/bin:$PATH
export LD_LIBRARY_PATH=/opt/pgi/linux86-64/16.5/lib:$MPI_HOME/lib:$CUDA_HOME/lib64:$LD_LIBRARY_PATH

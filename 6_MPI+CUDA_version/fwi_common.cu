#include "fwi_common.h"

#include <hwloc.h>
#include <hwloc/cudart.h>

#define CUDA_CHECK(call) { gpu_assert((call), __FILE__, __LINE__); }
inline void gpu_assert(cudaError_t code, const char* file, int line)
{
    if (code != cudaSuccess)
    {
        fprintf(stderr, "CUDA ERROR: %s:%d, ", file, line);
        fprintf(stderr, "code: %d, reason: %s\n", code, cudaGetErrorString(code));
        exit(code);
    }
}

extern "C"
int select_gpu_and_pin_proc(int rank, int local_rank)
{
    hwloc_topology_t topo;
    // Load hardware topology of all PCI devices in this node
    hwloc_topology_init(&topo);
    hwloc_topology_set_flags(topo, HWLOC_TOPOLOGY_FLAG_IO_DEVICES);
    hwloc_topology_load(topo);
    
    // choose a GPU based on MPI local rank
    int ngpus;
    CUDA_CHECK( cudaGetDeviceCount(&ngpus) );
    printf("acc_get_num_devices found %d devices (of type 'acc_device_nvidia')\n", ngpus);
    int device = local_rank % ngpus;
    //acc_set_gpu(device);
    cudaSetDevice(device);

    // Iterate through all CPU cores that are physically close to the selected GPU.
    // evenly distributing all processes across cores using local_rank
    hwloc_cpuset_t cpuset = hwloc_bitmap_alloc();
    hwloc_cudart_get_device_cpuset(topo, device, cpuset);
    int match = 0;
    int cpu;
    int i;
    hwloc_bitmap_foreach_begin(i, cpuset);
        if (match == local_rank)
        {
            cpu = i;
            break;
        }
        ++match;
    hwloc_bitmap_foreach_end();
    
    // Bind this process to the selected CPU
    hwloc_cpuset_t selcpu = hwloc_bitmap_alloc();
    hwloc_bitmap_set(selcpu, cpu);
    hwloc_set_cpubind(topo, selcpu, 0);
    
    // cleanup
    hwloc_bitmap_free(selcpu);
    hwloc_bitmap_free(cpuset);
    hwloc_topology_destroy(topo);
    
    char hostname[256];
    gethostname( hostname, sizeof(hostname) );
    cpu = sched_getcpu();
    printf("MPI rank %d [local rank %d] using GPU %d and CPU %d on host %s\n",
        rank, local_rank, device, cpu, hostname);
    return 0;
}

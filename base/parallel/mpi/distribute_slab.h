#pragma once

#include "base/field/field2.h"

#include <mpi.h>

namespace MPIUtils
{
    /**
     * Redistribute a matrix from one communicator's decomposition to another.
     * Assumes slab decomposition along the x-axis for both comm_src and comm_dest.
     * Since field2 is row-major (i*ny + j), x-slabs are contiguous in memory.
     */
    void redistribute_2d_slab_sync(field2* matrix_src, field2* matrix_dest, MPI_Comm comm_src, MPI_Comm comm_dest);
} // namespace MPIUtils
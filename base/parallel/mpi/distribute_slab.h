#pragma once

#include "base/field/field2.h"

#include <mpi.h>

namespace MPIUtils
{
    /**
     * @brief Redistribute a matrix from one communicator's decomposition to another.
     * Assumes slab decomposition along the x-axis for both comm_src and comm_dest.
     * Since field2 is row-major (i*ny + j), x-slabs are contiguous in memory.
     *
     * @param matrix_src is valid in comm_src, is nullptr in other ranks
     * @param matrix_dest is valid in comm_dest, is nullptr in other ranks
     * @param comm_src is valid in comm_src, is MPI_COMM_NULL in other ranks
     * @param comm_dest is valid in comm_dest, is MPI_COMM_NULL in other ranks
     */
    void redistribute_2d_slab_sync(field2* matrix_src, field2* matrix_dest, MPI_Comm comm_src, MPI_Comm comm_dest);
} // namespace MPIUtils
#pragma once

#include "base/field/field2.h"

#include <mpi.h>

namespace MPIUtils
{
    /**
     * @brief Redistributes a 1D double buffer from one communicator's decomposition to another.
     * Uses MPI_COMM_WORLD as the transport layer to bridge different sub-communicators.
     * * @param buffer_src  is valid in comm_src, is nullptr in other ranks.
     * @param buffer_dest is valid in comm_dest, is nullptr in other ranks.
     * @param n_src       Local length of the buffer in comm_src.
     * @param n_dest      Local length of the buffer in comm_dest.
     * @param comm_src is valid in comm_src, is MPI_COMM_NULL in other ranks
     * @param comm_dest is valid in comm_dest, is MPI_COMM_NULL in other ranks
     */
    void redistribute_slab(double*  buffer_src,
                           double*  buffer_dest,
                           int      n_src,
                           int      n_dest,
                           MPI_Comm comm_src,
                           MPI_Comm comm_dest);

    /**
     * @brief Redistribute a matrix from one communicator's decomposition to another.
     * Assumes slab decomposition along the x-axis for both comm_src and comm_dest.
     * Since field2 is row-major (i*ny + j), x-slabs are contiguous in memory.
     * matrix_src is [nx_slab_src, ny], matrix_dest is [nx_slab_dest, ny].
     * Sum of nx_slab_src in comm_src = sum of nx_slab_dest in comm_dest.
     *
     * @param matrix_src is valid in comm_src, is nullptr in other ranks.
     * @param matrix_dest is valid in comm_dest, is nullptr in other ranks.
     * @param comm_src is valid in comm_src, is MPI_COMM_NULL in other ranks
     * @param comm_dest is valid in comm_dest, is MPI_COMM_NULL in other ranks
     */
    void redistribute_slab(field2* matrix_src, field2* matrix_dest, MPI_Comm comm_src, MPI_Comm comm_dest);
} // namespace MPIUtils
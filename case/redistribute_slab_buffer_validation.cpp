#include "base/parallel/mpi/distribute_slab.h"
#include "base/parallel/mpi/mpi_misc.h"
#include <iomanip>
#include <iostream>
#include <vector>

int main(int argc, char** argv)
{
    MPI_Init(&argc, &argv);
    int world_rank, world_size;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    if (world_size < 4)
    {
        if (world_rank == 0)
            std::cerr << "This test requires at least 4 processes." << std::endl;
        MPI_Finalize();
        return 1;
    }

    const int n_total = 20; // Total length of the 1D buffer

    // --- Communicator Setup ---
    // Source: Ranks 0, 1, 2
    int      color_src = (world_rank <= 2) ? 0 : MPI_UNDEFINED;
    MPI_Comm comm_src;
    MPI_Comm_split(MPI_COMM_WORLD, color_src, world_rank, &comm_src);

    // Destination: Ranks 2, 3
    int      color_dest = (world_rank >= 2) ? 0 : MPI_UNDEFINED;
    MPI_Comm comm_dest;
    MPI_Comm_split(MPI_COMM_WORLD, color_dest, world_rank, &comm_dest);

    // --- Source Allocation & Initialization ---
    double* buffer_src = nullptr;
    int     n_src      = 0;
    if (comm_src != MPI_COMM_NULL)
    {
        int r, s;
        MPI_Comm_rank(comm_src, &r);
        MPI_Comm_size(comm_src, &s);

        n_src         = MPIUtils::get_slab_length(n_total, r, s);
        int start_idx = MPIUtils::get_slab_displacement(n_total, r, s);

        buffer_src = new double[n_src];
        for (int i = 0; i < n_src; ++i)
        {
            // Value is simply the global index for easy verification
            buffer_src[i] = static_cast<double>(start_idx + i);
        }
    }

    // --- Destination Allocation ---
    double* buffer_dest = nullptr;
    int     n_dest      = 0;
    if (comm_dest != MPI_COMM_NULL)
    {
        int r, s;
        MPI_Comm_rank(comm_dest, &r);
        MPI_Comm_size(comm_dest, &s);

        n_dest = MPIUtils::get_slab_length(n_total, r, s);

        buffer_dest = new double[n_dest];
        // Fill with dummy value to verify overwrite
        for (int i = 0; i < n_dest; ++i)
            buffer_dest[i] = -1.0;
    }

    // Ensure all allocations are done before proceeding
    MPI_Barrier(MPI_COMM_WORLD);

    // --- CALL REDISTRIBUTE ---
    // Rank 0, 1: Valid src, null dest
    // Rank 2:    Valid src, valid dest (Overlap point)
    // Rank 3:    Null src, valid dest
    MPIUtils::redistribute_slab(buffer_src, buffer_dest, n_src, n_dest, comm_src, comm_dest);

    // --- Verification ---
    if (comm_dest != MPI_COMM_NULL)
    {
        int r_dest;
        MPI_Comm_rank(comm_dest, &r_dest);
        int start_idx = MPIUtils::get_slab_displacement(n_total, r_dest, 2); // size_dest is 2

        std::cout << "[Rank " << world_rank << " / Dest Rank " << r_dest << "] received: ";
        for (int i = 0; i < n_dest; ++i)
        {
            std::cout << buffer_dest[i] << (i == n_dest - 1 ? "" : ", ");
        }
        std::cout << std::endl;
    }

    // --- Cleanup ---
    if (buffer_src)
        delete[] buffer_src;
    if (buffer_dest)
        delete[] buffer_dest;

    if (comm_src != MPI_COMM_NULL)
        MPI_Comm_free(&comm_src);
    if (comm_dest != MPI_COMM_NULL)
        MPI_Comm_free(&comm_dest);

    MPI_Finalize();
    return 0;
}
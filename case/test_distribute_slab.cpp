#include "base/parallel/mpi/distribute_slab.h"
#include "base/parallel/mpi/mpi_misc.h"
#include <iostream>

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

    const int nx_full = 12;
    const int ny_full = 4;

    // Create comm_src (Ranks 0, 1, 2)
    int      color_src = (world_rank <= 2) ? 0 : MPI_UNDEFINED;
    MPI_Comm comm_src;
    MPI_Comm_split(MPI_COMM_WORLD, color_src, world_rank, &comm_src);

    // Create comm_dest (Ranks 2, 3)
    int      color_dest = (world_rank >= 2) ? 0 : MPI_UNDEFINED;
    MPI_Comm comm_dest;
    MPI_Comm_split(MPI_COMM_WORLD, color_dest, world_rank, &comm_dest);

    // --- Source Allocation ---
    field2* ptr_src = nullptr;
    if (comm_src != MPI_COMM_NULL)
    {
        ptr_src = new field2();
        int r, s;
        MPI_Comm_rank(comm_src, &r);
        MPI_Comm_size(comm_src, &s);
        int local_nx = MPIUtils::get_slab_length(nx_full, r, s);
        int start_nx = MPIUtils::get_slab_displacement(nx_full, r, s);
        ptr_src->init(local_nx, ny_full, "Source_Field");
        for (int i = 0; i < local_nx; ++i)
            for (int j = 0; j < ny_full; ++j)
                (*ptr_src)(i, j) = 100.0 * (start_nx + i) + j;
    }

    // --- Destination Allocation ---
    field2* ptr_dest = nullptr;
    if (comm_dest != MPI_COMM_NULL)
    {
        ptr_dest = new field2();
        int r, s;
        MPI_Comm_rank(comm_dest, &r);
        MPI_Comm_size(comm_dest, &s);
        int local_nx = MPIUtils::get_slab_length(nx_full, r, s);
        ptr_dest->init(local_nx, ny_full, "Dest_Field");
        ptr_dest->clear(-1.0);
    }

    MPI_Barrier(MPI_COMM_WORLD);

    // --- CALL WITH POINTERS ---
    // Process 0,1: (field, nullptr)
    // Process 2:   (field, field)
    // Process 3:   (nullptr, field)
    MPIUtils::redistribute_2d_slab_sync(ptr_src, ptr_dest, comm_src, comm_dest);

    // --- Verification ---
    if (ptr_dest)
    {
        // ... (printing logic same as your original, just use ptr_dest->)
        std::cout << "Rank " << world_rank << " received data." << std::endl;
        ptr_dest->print();
    }

    // Cleanup
    delete ptr_src;
    delete ptr_dest;
    if (comm_src != MPI_COMM_NULL)
        MPI_Comm_free(&comm_src);
    if (comm_dest != MPI_COMM_NULL)
        MPI_Comm_free(&comm_dest);

    MPI_Finalize();
    return 0;
}
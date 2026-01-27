#include "base/parallel/mpi/mpi_misc.h"
#include "base/parallel/mpi/transpose_slab.h"

int main(int argc, char* argv[])
{
    MPI_Init(&argc, &argv);

    int nx = 4;
    int ny = 3;

    int mpi_rank = 0, mpi_size = 1;

    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);

    int nx_slab = MPIUtils::get_slab_length(nx, mpi_rank, mpi_size);
    int nx_disp = MPIUtils::get_slab_displacement(nx, mpi_rank, mpi_size);
    int ny_slab = MPIUtils::get_slab_length(ny, mpi_rank, mpi_size);

    field2 f(nx_slab, ny), f_T(ny_slab, nx);
    for (int i = 0; i < nx_slab; i++)
    {
        for (int j = 0; j < ny; j++)
        {
            f(i, j) = (i + nx_disp) * ny + j + 1;
        }
    }

    for (int i = 0; i < mpi_size; i++)
    {
        if (i == mpi_rank)
        {
            std::cout << "rank " << mpi_rank << std::endl;
            f.print();
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }

    if (mpi_rank == 0)
        std::cout << "----------------------------" << std::endl;

    double* buf_slab_x = new double[nx_slab * ny];
    double* buf_slab_y = new double[ny_slab * nx];

    MPIUtils::transpose_2d_slab_sync(f, buf_slab_x, f_T, buf_slab_y, MPI_COMM_WORLD);

    for (int i = 0; i < mpi_size; i++)
    {
        if (i == mpi_rank)
        {
            std::cout << "rank " << mpi_rank << std::endl;
            f_T.print();
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }

    MPI_Finalize();

    return 0;
}
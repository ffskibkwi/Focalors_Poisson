#include "transpose_slab.h"

namespace MPIUtils
{
    void transpose_2d_slab_sync(field2&  matrix,
                                double*  send_buffer,
                                field2&  matrix_T,
                                double*  recv_buffer,
                                MPI_Comm communicator)
    {
        // scatter along y
        // every part is [nx, ny * p / mpi_size ~ ny * (p + 1) / mpi_size]

        int mpi_size = 1;
        MPI_Comm_size(communicator, &mpi_size);

        const int nx_block = matrix.get_nx();
        const int ny       = matrix.get_ny();
        const int nx       = nx_block * mpi_size;
        const int ny_block = ny / mpi_size;

        OPENMP_PARALLEL_FOR()
        for (int p = 0; p < mpi_size; ++p)
        {
            for (int i = 0; i < nx_block; ++i)
            {
                for (int j = p * ny_block; j < (p + 1) * ny_block; ++j)
                {
                    send_buffer[p * nx_block * ny_block + i * ny_block + (j - p * ny_block)] = matrix(i, j);
                }
            }
        }

        MPI_Alltoall(
            send_buffer, nx_block * ny_block, MPI_DOUBLE, recv_buffer, nx_block * ny_block, MPI_DOUBLE, communicator);

        int recv_index = 0;
        for (int p = 0; p < mpi_size; ++p)
        {
            for (int i = p * nx_block; i < (p + 1) * nx_block; ++i)
            {
                for (int j = 0; j < ny_block; ++j)
                {
                    matrix_T(j, i) = recv_buffer[recv_index++];
                }
            }
        }
    }

    void transpose_2d_slab_async_begin(field2&      matrix,
                                       double*      send_buffer,
                                       double*      recv_buffer,
                                       MPI_Comm     communicator,
                                       MPI_Request* request)
    {
        // scatter along y
        // every part is [nx, ny * p / mpi_size ~ ny * (p + 1) / mpi_size]

        int mpi_size = 1;
        MPI_Comm_size(communicator, &mpi_size);

        const int nx_block = matrix.get_nx();
        const int ny       = matrix.get_ny();
        const int nx       = nx_block * mpi_size;
        const int ny_block = ny / mpi_size;

        OPENMP_PARALLEL_FOR()
        for (int p = 0; p < mpi_size; ++p)
        {
            for (int i = 0; i < nx_block; ++i)
            {
                for (int j = p * ny_block; j < (p + 1) * ny_block; ++j)
                {
                    send_buffer[p * nx_block * ny_block + i * ny_block + (j - p * ny_block)] = matrix(i, j);
                }
            }
        }

        MPI_Ialltoall(send_buffer,
                      nx_block * ny_block,
                      MPI_DOUBLE,
                      recv_buffer,
                      nx_block * ny_block,
                      MPI_DOUBLE,
                      communicator,
                      request);
    }

    void transpose_2d_slab_async_end(field2& matrix_T, double* recv_buffer, MPI_Comm communicator, MPI_Request* request)
    {
        int mpi_size = 1;
        MPI_Comm_size(communicator, &mpi_size);

        const int ny_block = matrix_T.get_nx();
        const int nx       = matrix_T.get_ny();
        const int ny       = ny_block * mpi_size;
        const int nx_block = nx / mpi_size;

        MPI_Wait(request, MPI_STATUS_IGNORE);

        int recv_index = 0;
        for (int p = 0; p < mpi_size; ++p)
        {
            for (int i = p * nx_block; i < (p + 1) * nx_block; ++i)
            {
                for (int j = 0; j < ny_block; ++j)
                {
                    matrix_T(j, i) = recv_buffer[recv_index++];
                }
            }
        }
    }

    void transpose_3d_slab_shift_right_sync(field3&  matrix,
                                            double*  send_buffer,
                                            field3&  matrix_T,
                                            double*  recv_buffer,
                                            MPI_Comm communicator)
    {
        int mpi_size = 1;
        MPI_Comm_size(communicator, &mpi_size);

        const int nx_block = matrix.get_nx();
        const int ny       = matrix.get_ny();
        const int nz       = matrix.get_nz();

        const int nz_block = nz / mpi_size;
        const int nx       = nx_block * mpi_size;

        OPENMP_PARALLEL_FOR()
        for (int i = 0; i < nx_block; ++i)
        {
            for (int j = 0; j < ny; ++j)
            {
                for (int k = 0; k < nz; ++k)
                {
                    // (nx_block, ny, nz) -> (nz, nx_block, ny)
                    send_buffer[k * nx_block * ny + i * ny + j] = matrix(i, j, k);
                }
            }
        }

        MPI_Alltoall(send_buffer,
                     nz_block * nx_block * ny,
                     MPI_DOUBLE,
                     recv_buffer,
                     nz_block * nx_block * ny,
                     MPI_DOUBLE,
                     communicator);

        OPENMP_PARALLEL_FOR()
        for (int p = 0; p < mpi_size; ++p)
        {
            for (int i = 0; i < nz_block; ++i)
            {
                for (int j = p * nx_block; j < (p + 1) * nx_block; ++j)
                {
                    for (int k = 0; k < ny; ++k)
                    {
                        matrix_T(i, j, k) =
                            recv_buffer[p * nz_block * nx_block * ny + i * nx_block * ny + (j - p * nx_block) * ny + k];
                    }
                }
            }
        }
    }

    void transpose_3d_slab_shift_right_async_begin(field3&      matrix,
                                                   double*      send_buffer,
                                                   double*      recv_buffer,
                                                   MPI_Comm     communicator,
                                                   MPI_Request* request)
    {
        int mpi_size = 1;
        MPI_Comm_size(communicator, &mpi_size);

        const int nx_block = matrix.get_nx();
        const int ny       = matrix.get_ny();
        const int nz       = matrix.get_nz();

        const int nz_block = nz / mpi_size;
        const int nx       = nx_block * mpi_size;

        OPENMP_PARALLEL_FOR()
        for (int i = 0; i < nx_block; ++i)
        {
            for (int j = 0; j < ny; ++j)
            {
                for (int k = 0; k < nz; ++k)
                {
                    // (nx_block, ny, nz) -> (nz, nx_block, ny)
                    send_buffer[k * nx_block * ny + i * ny + j] = matrix(i, j, k);
                }
            }
        }

        MPI_Ialltoall(send_buffer,
                      nz_block * nx_block * ny,
                      MPI_DOUBLE,
                      recv_buffer,
                      nz_block * nx_block * ny,
                      MPI_DOUBLE,
                      communicator,
                      request);
    }

    void transpose_3d_slab_shift_right_async_end(field3&      matrix_T,
                                                 double*      recv_buffer,
                                                 MPI_Comm     communicator,
                                                 MPI_Request* request)
    {
        int mpi_size = 1;
        MPI_Comm_size(communicator, &mpi_size);

        const int nz_block = matrix_T.get_nx();
        const int nx       = matrix_T.get_ny();
        const int ny       = matrix_T.get_nz();

        const int nx_block = nx / mpi_size;
        const int nz       = nz_block * mpi_size;

        MPI_Wait(request, MPI_STATUS_IGNORE);

        OPENMP_PARALLEL_FOR()
        for (int p = 0; p < mpi_size; ++p)
        {
            for (int i = 0; i < nz_block; ++i)
            {
                for (int j = p * nx_block; j < (p + 1) * nx_block; ++j)
                {
                    for (int k = 0; k < ny; ++k)
                    {
                        matrix_T(i, j, k) =
                            recv_buffer[p * nz_block * nx_block * ny + i * nx_block * ny + (j - p * nx_block) * ny + k];
                    }
                }
            }
        }
    }

    void transpose_3d_slab_shift_left_sync(field3&  matrix,
                                           double*  send_buffer,
                                           field3&  matrix_T,
                                           double*  recv_buffer,
                                           MPI_Comm communicator)
    {
        int mpi_size = 1;
        MPI_Comm_size(communicator, &mpi_size);

        const int nx_block = matrix.get_nx();
        const int ny       = matrix.get_ny();
        const int nz       = matrix.get_nz();

        const int ny_block = ny / mpi_size;
        const int nx       = nx_block * mpi_size;

        OPENMP_PARALLEL_FOR()
        for (int i = 0; i < nx_block; ++i)
        {
            for (int j = 0; j < ny; ++j)
            {
                for (int k = 0; k < nz; ++k)
                {
                    // (nx_block, ny, nz) -> (ny, nz, nx_block)
                    send_buffer[j * nz * nx_block + k * nx_block + i] = matrix(i, j, k);
                }
            }
        }

        MPI_Alltoall(send_buffer,
                     ny_block * nz * nx_block,
                     MPI_DOUBLE,
                     recv_buffer,
                     ny_block * nz * nx_block,
                     MPI_DOUBLE,
                     communicator);

        OPENMP_PARALLEL_FOR()
        for (int p = 0; p < mpi_size; ++p)
        {
            for (int i = 0; i < ny_block; ++i)
            {
                for (int j = 0; j < nz; ++j)
                {
                    for (int k = p * nx_block; k < (p + 1) * nx_block; ++k)
                    {
                        matrix_T(i, j, k) = recv_buffer[p * ny_block * nz * nx_block + i * nz * nx_block +
                                                        j * nx_block + (k - p * nx_block)];
                    }
                }
            }
        }
    }

    void transpose_3d_slab_shift_left_async_begin(field3&      matrix,
                                                  double*      send_buffer,
                                                  double*      recv_buffer,
                                                  MPI_Comm     communicator,
                                                  MPI_Request* request)
    {
        int mpi_size = 1;
        MPI_Comm_size(communicator, &mpi_size);

        const int nx_block = matrix.get_nx();
        const int ny       = matrix.get_ny();
        const int nz       = matrix.get_nz();

        const int ny_block = ny / mpi_size;
        const int nx       = nx_block * mpi_size;

        OPENMP_PARALLEL_FOR()
        for (int i = 0; i < nx_block; ++i)
        {
            for (int j = 0; j < ny; ++j)
            {
                for (int k = 0; k < nz; ++k)
                {
                    // (nx_block, ny, nz) -> (ny, nz, nx_block)
                    send_buffer[j * nz * nx_block + k * nx_block + i] = matrix(i, j, k);
                }
            }
        }

        MPI_Ialltoall(send_buffer,
                      ny_block * nz * nx_block,
                      MPI_DOUBLE,
                      recv_buffer,
                      ny_block * nz * nx_block,
                      MPI_DOUBLE,
                      communicator,
                      request);
    }

    void transpose_3d_slab_shift_left_async_end(field3&      matrix_T,
                                                double*      recv_buffer,
                                                MPI_Comm     communicator,
                                                MPI_Request* request)
    {
        int mpi_size = 1;
        MPI_Comm_size(communicator, &mpi_size);

        const int ny_block = matrix_T.get_nx();
        const int nz       = matrix_T.get_ny();
        const int nx       = matrix_T.get_nz();

        const int nx_block = nx / mpi_size;
        const int ny       = ny_block * mpi_size;

        MPI_Wait(request, MPI_STATUS_IGNORE);

        OPENMP_PARALLEL_FOR()
        for (int p = 0; p < mpi_size; ++p)
        {
            for (int i = 0; i < ny_block; ++i)
            {
                for (int j = 0; j < nz; ++j)
                {
                    for (int k = p * nx_block; k < (p + 1) * nx_block; ++k)
                    {
                        matrix_T(i, j, k) = recv_buffer[p * ny_block * nz * nx_block + i * nz * nx_block +
                                                        j * nx_block + (k - p * nx_block)];
                    }
                }
            }
        }
    }
} // namespace MPIUtils
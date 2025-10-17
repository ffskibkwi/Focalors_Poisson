#include "transpose_pencil.h"

namespace MPIUtils
{
    void transpose_3d_pencil_shift_right_sync(field3&  matrix,
                                              double*  send_buffer,
                                              field3&  matrix_T,
                                              double*  recv_buffer,
                                              MPI_Comm communicator)
    {
        int mpi_size = 1;
        MPI_Comm_size(communicator, &mpi_size);

        const int nx_block0 = matrix.get_nx();
        const int ny_block1 = matrix.get_ny();
        const int nz        = matrix.get_nz();

        const int nz_block1 = nz / mpi_size;
        const int ny        = ny_block1 * mpi_size;

        OPENMP_PARALLEL_FOR()
        for (int i = 0; i < nx_block0; ++i)
        {
            for (int j = 0; j < ny_block1; ++j)
            {
                for (int k = 0; k < nz; ++k)
                {
                    // (nx_block0, ny_block1, nz) -> (nz, nx_block0, ny_block1)
                    send_buffer[k * nx_block0 * ny_block1 + i * ny_block1 + j] = matrix(i, j, k);
                }
            }
        }

        MPI_Alltoall(send_buffer,
                     nz_block1 * nx_block0 * ny_block1,
                     MPI_DOUBLE,
                     recv_buffer,
                     nz_block1 * nx_block0 * ny_block1,
                     MPI_DOUBLE,
                     communicator);

        OPENMP_PARALLEL_FOR()
        for (int p = 0; p < mpi_size; ++p)
        {
            for (int i = 0; i < nz_block1; ++i)
            {
                for (int j = 0; j < nx_block0; ++j)
                {
                    for (int k = p * ny_block1; k < (p + 1) * ny_block1; ++k)
                    {
                        matrix_T(i, j, k) =
                            recv_buffer[p * nz_block1 * nx_block0 * ny_block1 + i * nx_block0 * ny_block1 +
                                        j * ny_block1 + (k - p * ny_block1)];
                    }
                }
            }
        }
    }

    void transpose_3d_pencil_shift_right_async_begin(field3&      matrix,
                                                     double*      send_buffer,
                                                     double*      recv_buffer,
                                                     MPI_Comm     communicator,
                                                     MPI_Request* request)
    {
        int mpi_size = 1;
        MPI_Comm_size(communicator, &mpi_size);

        const int nx_block0 = matrix.get_nx();
        const int ny_block1 = matrix.get_ny();
        const int nz        = matrix.get_nz();

        const int nz_block1 = nz / mpi_size;
        const int ny        = ny_block1 * mpi_size;

        OPENMP_PARALLEL_FOR()
        for (int i = 0; i < nx_block0; ++i)
        {
            for (int j = 0; j < ny_block1; ++j)
            {
                for (int k = 0; k < nz; ++k)
                {
                    // (nx_block0, ny_block1, nz) -> (nz, nx_block0, ny_block1)
                    send_buffer[k * nx_block0 * ny_block1 + i * ny_block1 + j] = matrix(i, j, k);
                }
            }
        }

        MPI_Ialltoall(send_buffer,
                      nz_block1 * nx_block0 * ny_block1,
                      MPI_DOUBLE,
                      recv_buffer,
                      nz_block1 * nx_block0 * ny_block1,
                      MPI_DOUBLE,
                      communicator,
                      request);
    }

    void transpose_3d_pencil_shift_right_async_end(field3&      matrix_T,
                                                   double*      recv_buffer,
                                                   MPI_Comm     communicator,
                                                   MPI_Request* request)
    {
        int mpi_size = 1;
        MPI_Comm_size(communicator, &mpi_size);

        const int nz_block1 = matrix_T.get_nx();
        const int nx_block0 = matrix_T.get_ny();
        const int ny        = matrix_T.get_nz();

        const int ny_block1 = ny / mpi_size;
        const int nz        = nz_block1 * mpi_size;

        MPI_Wait(request, MPI_STATUS_IGNORE);

        OPENMP_PARALLEL_FOR()
        for (int p = 0; p < mpi_size; ++p)
        {
            for (int i = 0; i < nz_block1; ++i)
            {
                for (int j = 0; j < nx_block0; ++j)
                {
                    for (int k = p * ny_block1; k < (p + 1) * ny_block1; ++k)
                    {
                        matrix_T(i, j, k) =
                            recv_buffer[p * nz_block1 * nx_block0 * ny_block1 + i * nx_block0 * ny_block1 +
                                        j * ny_block1 + (k - p * ny_block1)];
                    }
                }
            }
        }
    }

    void transpose_3d_pencil_shift_left_sync(field3&  matrix,
                                             double*  send_buffer,
                                             field3&  matrix_T,
                                             double*  recv_buffer,
                                             MPI_Comm communicator)
    {
        int mpi_size = 1;
        MPI_Comm_size(communicator, &mpi_size);

        const int nx_block0 = matrix.get_nx();
        const int ny_block1 = matrix.get_ny();
        const int nz        = matrix.get_nz();

        const int nz_block0 = nz / mpi_size;
        const int nx        = nx_block0 * mpi_size;

        OPENMP_PARALLEL_FOR()
        for (int i = 0; i < nx_block0; ++i)
        {
            for (int j = 0; j < ny_block1; ++j)
            {
                for (int k = 0; k < nz; ++k)
                {
                    // (nx_block0, ny_block1, nz) -> (nz, ny_block1, nx_block0)
                    send_buffer[k * ny_block1 * nx_block0 + j * nx_block0 + i] = matrix(i, j, k);
                }
            }
        }

        MPI_Alltoall(send_buffer,
                     nz_block0 * ny_block1 * nx_block0,
                     MPI_DOUBLE,
                     recv_buffer,
                     nz_block0 * ny_block1 * nx_block0,
                     MPI_DOUBLE,
                     communicator);

        OPENMP_PARALLEL_FOR()
        for (int p = 0; p < mpi_size; ++p)
        {
            for (int i = 0; i < nz_block0; ++i)
            {
                for (int j = 0; j < ny_block1; ++j)
                {
                    for (int k = p * nx_block0; k < (p + 1) * nx_block0; ++k)
                    {
                        matrix_T(j, i, k) =
                            recv_buffer[p * nz_block0 * ny_block1 * nx_block0 + i * ny_block1 * nx_block0 +
                                        j * nx_block0 + (k - p * nx_block0)];
                    }
                }
            }
        }
    }

    void transpose_3d_pencil_shift_left_async_begin(field3&      matrix,
                                                    double*      send_buffer,
                                                    double*      recv_buffer,
                                                    MPI_Comm     communicator,
                                                    MPI_Request* request)
    {
        int mpi_size = 1;
        MPI_Comm_size(communicator, &mpi_size);

        const int nx_block0 = matrix.get_nx();
        const int ny_block1 = matrix.get_ny();
        const int nz        = matrix.get_nz();

        const int nz_block0 = nz / mpi_size;
        const int nx        = nx_block0 * mpi_size;

        OPENMP_PARALLEL_FOR()
        for (int i = 0; i < nx_block0; ++i)
        {
            for (int j = 0; j < ny_block1; ++j)
            {
                for (int k = 0; k < nz; ++k)
                {
                    // (nx_block0, ny_block1, nz) -> (nz, ny_block1, nx_block0)
                    send_buffer[k * ny_block1 * nx_block0 + j * nx_block0 + i] = matrix(i, j, k);
                }
            }
        }

        MPI_Ialltoall(send_buffer,
                      nz_block0 * ny_block1 * nx_block0,
                      MPI_DOUBLE,
                      recv_buffer,
                      nz_block0 * ny_block1 * nx_block0,
                      MPI_DOUBLE,
                      communicator,
                      request);
    }

    void transpose_3d_pencil_shift_left_async_end(field3&      matrix_T,
                                                  double*      recv_buffer,
                                                  MPI_Comm     communicator,
                                                  MPI_Request* request)
    {
        int mpi_size = 1;
        MPI_Comm_size(communicator, &mpi_size);

        const int ny_block1 = matrix_T.get_nx();
        const int nz_block0 = matrix_T.get_ny();
        const int nx        = matrix_T.get_nz();

        const int nx_block0 = nx / mpi_size;
        const int nz        = nz_block0 * mpi_size;

        MPI_Wait(request, MPI_STATUS_IGNORE);

        OPENMP_PARALLEL_FOR()
        for (int p = 0; p < mpi_size; ++p)
        {
            for (int i = 0; i < nz_block0; ++i)
            {
                for (int j = 0; j < ny_block1; ++j)
                {
                    for (int k = p * nx_block0; k < (p + 1) * nx_block0; ++k)
                    {
                        matrix_T(j, i, k) =
                            recv_buffer[p * nz_block0 * ny_block1 * nx_block0 + i * ny_block1 * nx_block0 +
                                        j * nx_block0 + (k - p * nx_block0)];
                    }
                }
            }
        }
    }
} // namespace MPIUtils
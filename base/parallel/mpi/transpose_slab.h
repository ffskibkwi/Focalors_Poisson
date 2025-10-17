#pragma once

#include "pch.h"

#include <mpi.h>

/**
 * @brief To perform global redistributions, there is two steps necessary(not necessarily in this order):
 *
 * 1. Contiguous data communication using collective all-to-all operations
 *
 * 2. Local data rearrangements or transpose operations, also referred as remappings
 */
namespace MPIUtils
{
    /**
     * @brief Transpose a distributed 2d matrix from xy to yx.
     *
     * @param matrix [nx/P, ny].
     * @param send_buffer [nx/P * ny].
     * @param matrix_T [ny/P, nx].
     * @param recv_buffer [nx/P * ny].
     * @param communicator The communicator in which the standard send takes place.
     */
    void transpose_2d_slab_sync(field2&  matrix,
                                double*  send_buffer,
                                field2&  matrix_T,
                                double*  recv_buffer,
                                MPI_Comm communicator);

    /**
     * @brief Transpose a distributed 2d matrix from xy to yx.
     *
     * @param matrix [nx/P, ny].
     * @param send_buffer [nx/P * ny].
     * @param recv_buffer [nx/P * ny].
     * @param communicator The communicator in which the standard send takes place.
     * @param request The variable in which store the handler on the non-blocking operation.
     */
    void transpose_2d_slab_async_begin(field2&      matrix,
                                       double*      send_buffer,
                                       double*      recv_buffer,
                                       MPI_Comm     communicator,
                                       MPI_Request* request);

    /**
     * @brief Transpose a distributed 2d matrix from xy to yx.
     *
     * @param matrix_T [ny/P, nx].
     * @param recv_buffer [nx/P * ny].
     * @param communicator The communicator in which the standard send takes place.
     * @param request The variable in which store the handler on the non-blocking operation.
     */
    void
    transpose_2d_slab_async_end(field2& matrix_T, double* recv_buffer, MPI_Comm communicator, MPI_Request* request);

    /**
     * @brief Transpose a distributed 3d matrix from xyz to zxy.
     *
     * @param matrix [nx/P, ny, nz].
     * @param send_buffer [nx/P * ny * nz].
     * @param matrix_T [nz/P, nx, ny].
     * @param recv_buffer [nx/P * ny * nz].
     * @param communicator The communicator in which the standard send takes place.
     */
    void transpose_3d_slab_shift_right_sync(field3&  matrix,
                                            double*  send_buffer,
                                            field3&  matrix_T,
                                            double*  recv_buffer,
                                            MPI_Comm communicator);

    /**
     * @brief Transpose a distributed 3d matrix from xyz to zxy.
     *
     * @param matrix [nx/P, ny, nz].
     * @param send_buffer [nx/P * ny * nz].
     * @param recv_buffer [nx/P * ny * nz].
     * @param communicator The communicator in which the standard send takes place.
     * @param request The variable in which store the handler on the non-blocking operation.
     */
    void transpose_3d_slab_shift_right_async_begin(field3&      matrix,
                                                   double*      send_buffer,
                                                   double*      recv_buffer,
                                                   MPI_Comm     communicator,
                                                   MPI_Request* request);

    /**
     * @brief Transpose a distributed 3d matrix from xyz to zxy.
     *
     * @param matrix_T [nz/P, nx, ny].
     * @param recv_buffer [nx/P * ny * nz].
     * @param communicator The communicator in which the standard send takes place.
     * @param request The variable in which store the handler on the non-blocking operation.
     */
    void transpose_3d_slab_shift_right_async_end(field3&      matrix_T,
                                                 double*      recv_buffer,
                                                 MPI_Comm     communicator,
                                                 MPI_Request* request);

    /**
     * @brief Transpose a distributed 3d matrix from xyz to yzx.
     *
     * @param matrix [nx/P, ny, nz].
     * @param send_buffer [nx/P * ny * nz].
     * @param matrix_T [ny/P, nz, nx].
     * @param recv_buffer [nx/P * ny * nz].
     * @param communicator The communicator in which the standard send takes place.
     */
    void transpose_3d_slab_shift_left_sync(field3&  matrix,
                                           double*  send_buffer,
                                           field3&  matrix_T,
                                           double*  recv_buffer,
                                           MPI_Comm communicator);

    /**
     * @brief Transpose a distributed 3d matrix from xyz to yzx.
     *
     * @param matrix [nx/P, ny, nz].
     * @param send_buffer [nx/P * ny * nz].
     * @param recv_buffer [nx/P * ny * nz].
     * @param communicator The communicator in which the standard send takes place.
     * @param request The variable in which store the handler on the non-blocking operation.
     */
    void transpose_3d_slab_shift_left_async_begin(field3&      matrix,
                                                  double*      send_buffer,
                                                  double*      recv_buffer,
                                                  MPI_Comm     communicator,
                                                  MPI_Request* request);

    /**
     * @brief Transpose a distributed 3d matrix from xyz to yzx.
     *
     * @param matrix_T [ny/P, nz, nx].
     * @param recv_buffer [nx/P * ny * nz].
     * @param communicator The communicator in which the standard send takes place.
     * @param request The variable in which store the handler on the non-blocking operation.
     */
    void transpose_3d_slab_shift_left_async_end(field3&      matrix_T,
                                                double*      recv_buffer,
                                                MPI_Comm     communicator,
                                                MPI_Request* request);
} // namespace MPIUtils
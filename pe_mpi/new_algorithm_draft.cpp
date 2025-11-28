#include <mpi.h>
#include <vector>
#include <iostream>
#include <cmath>
#include <algorithm>
#include <cstdio>
#include <string>

// Helper to print distributed matrix
void print_distributed_matrix(const std::vector<double>& local_A, int local_rows, int M, int rank, int size, const std::string& label) {
    std::vector<int> all_rows(size);
    MPI_Gather(&local_rows, 1, MPI_INT, all_rows.data(), 1, MPI_INT, 0, MPI_COMM_WORLD);

    std::vector<double> global_A;
    std::vector<int> recvcounts;
    std::vector<int> displs;

    if (rank == 0) {
        recvcounts.resize(size);
        displs.resize(size);
        int offset = 0;
        for (int i = 0; i < size; ++i) {
            recvcounts[i] = all_rows[i] * M;
            displs[i] = offset;
            offset += recvcounts[i];
        }
        global_A.resize(offset);
    }

    MPI_Gatherv(local_A.data(), local_rows * M, MPI_DOUBLE,
                global_A.data(), recvcounts.data(), displs.data(), MPI_DOUBLE,
                0, MPI_COMM_WORLD);

    if (rank == 0) {
        std::cout << "\n[" << label << "]\n";
        int total_rows = global_A.size() / M;
        for (int i = 0; i < total_rows; ++i) {
            printf("Row %2d: ", i);
            for (int j = 0; j < M; ++j) {
                printf("%6.2f ", global_A[i * M + j]);
            }
            printf("\n");
        }
        std::cout << std::endl;
    }
}

// Abstract vector function as required
// v: pointer to the vector data (size n)
// n: size of the vector
// para: parameter for the calculation
void vector_func(double* v, int n, double para)
{
    // Example implementation: v = A_temp * v
    // Constructing a temporary matrix A_temp implicitly
    // For demonstration, let's just do a simple scaling and neighbor interaction
    // to simulate a matrix-vector multiplication.
    
    std::vector<double> v_in(v, v + n);
    
    for (int i = 0; i < n; ++i) {
        double val = 0.0;
        // Diagonal
        val += para * v_in[i];
        // Off-diagonal (n-1, 1) and (n-1, -1) simulation
        if (i > 0) val += v_in[i-1];
        if (i < n - 1) val += v_in[i+1];
        
        v[i] = val;
    }
}

// Helper to transpose a local block
// src: rows x cols
// dst: cols x rows
void local_transpose(const std::vector<double>& src, std::vector<double>& dst, int rows, int cols)
{
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            dst[j * rows + i] = src[i * cols + j];
        }
    }
}

int main(int argc, char** argv)
{
    MPI_Init(&argc, &argv);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    int N = 32;
    int M = 32;
    int block_size = 4; // block is the col block

    // 1. Determine local row range
    // Simple block distribution of rows
    int rows_per_proc = N / size;
    int remainder = N % size;
    int i_start = rank * rows_per_proc + std::min(rank, remainder);
    int i_end = i_start + rows_per_proc + (rank < remainder ? 1 : 0);
    int local_rows = i_end - i_start;

    // Exchange local_rows to build recvcounts for Reduce_scatter
    std::vector<int> all_local_rows(size);
    MPI_Allgather(&local_rows, 1, MPI_INT, all_local_rows.data(), 1, MPI_INT, MPI_COMM_WORLD);

    std::vector<int> recvcounts(size);
    for(int r=0; r<size; ++r) {
        recvcounts[r] = all_local_rows[r] * block_size;
    }

    // 2. Initialize Matrix A (Local storage)
    // Stored in row-major: A[local_row_idx][global_col_idx]
    // Flattened: A[local_row_idx * M + global_col_idx]
    std::vector<double> A(local_rows * M);
    for (int i = 0; i < local_rows; i++)
    {
        int global_i = i_start + i;
        for (int j = 0; j < M; j++)
        {
            A[i * M + j] = 1.1 * global_i + 2.3 * j;
        }
    }

    print_distributed_matrix(A, local_rows, M, rank, size, "Initial Matrix A");

    // 3. Pipeline Setup
    int num_blocks = (M + block_size - 1) / block_size;
    
    // Double buffering structures
    // We need buffers for the "Compute" phase and "Comm" phase.
    // Flow: 
    // A_block (local) -> Transpose -> Padded_Block (B x N) -> vector_func -> Transpose -> SendBuf (N x B) -> Reduce_scatter -> RecvBuf (local_rows x B) -> A_block
    
    // We need persistent buffers for the async communication.
    // send_buffers: Input to Reduce_scatter. Size: N * block_size.
    // recv_buffers: Output of Reduce_scatter. Size: local_rows * block_size.
    std::vector<std::vector<double>> send_buffers(2, std::vector<double>(N * block_size));
    std::vector<std::vector<double>> recv_buffers(2, std::vector<double>(local_rows * block_size));
    
    MPI_Request reqs[2] = {MPI_REQUEST_NULL, MPI_REQUEST_NULL};

    // Helper lambda to compute a block
    auto compute_block = [&](int ib, int buf_idx) {
        int current_block_width = std::min(block_size, M - ib * block_size);
        
        // 1. Extract and Transpose local block to temp buffer
        // We want Padded_Block to be (block_size x N) so we can iterate over rows (which are columns of A)
        // But A is (local_rows x M).
        // Let's extract local block: (local_rows x current_block_width)
        // And place it into a (block_size x N) buffer at the correct column offset (i_start).
        
        // Initialize the full N*block_size buffer with zeros
        std::fill(send_buffers[buf_idx].begin(), send_buffers[buf_idx].end(), 0.0);
        
        // Fill the local part
        // send_buffers is treated as (block_size x N) temporarily for vector_func
        // But wait, for Reduce_scatter we need (N x block_size) layout (row-major N, col-major B) 
        // OR (block_size x N) depending on how we scatter.
        // We decided: Reduce_scatter splits a 1D array.
        // Rank r gets `recvcounts[r]` elements.
        // `recvcounts[r] = local_rows[r] * block_size`.
        // This implies the data must be arranged such that the chunk for rank r is contiguous.
        // Chunk for rank r corresponds to rows [i_start, i_end) of the RESULT matrix.
        // The result matrix is (N x block_size).
        // So we need the data in Row-Major order of the (N x block_size) matrix.
        // i.e. Row 0 (all cols), Row 1 (all cols)...
        
        // So `send_buffers` should be prepared as (N x block_size).
        // `vector_func` operates on columns of A (length N).
        // So we have `block_size` columns.
        
        // Strategy:
        // 1. Create a temp buffer `work_cols` (block_size x N).
        // 2. Fill local part: `work_cols[c][i_start...i_end]` from A.
        // 3. Apply `vector_func` on each row of `work_cols`.
        // 4. Transpose `work_cols` (block_size x N) -> `send_buffers[buf_idx]` (N x block_size).
        
        std::vector<double> work_cols(block_size * N, 0.0);
        
        for (int c = 0; c < current_block_width; ++c) {
            int global_j = ib * block_size + c;
            for (int i = 0; i < local_rows; ++i) {
                // A is row-major: A[i][global_j]
                // work_cols is row-major (c, global_i): work_cols[c * N + (i_start + i)]
                work_cols[c * N + (i_start + i)] = A[i * M + global_j];
            }
        }
        
        // Compute
        for (int c = 0; c < current_block_width; ++c) {
            double para = 3.03 * c; // Dummy parameter
            vector_func(&work_cols[c * N], N, para);
        }
        
        // Transpose to send_buffers (N x block_size)
        // work_cols: (block_size, N) -> send_buffers: (N, block_size)
        local_transpose(work_cols, send_buffers[buf_idx], block_size, N);
    };

    // Helper lambda to unpack a block
    auto unpack_block = [&](int ib, int buf_idx) {
        int current_block_width = std::min(block_size, M - ib * block_size);
        // recv_buffers[buf_idx] contains (local_rows x block_size) data.
        // We need to put it back into A.
        
        for (int i = 0; i < local_rows; ++i) {
            for (int c = 0; c < current_block_width; ++c) {
                int global_j = ib * block_size + c;
                // recv_buffers is row-major (local_rows, block_size)
                A[i * M + global_j] = recv_buffers[buf_idx][i * block_size + c];
            }
        }
    };

    // 4. Execution Loop
    // Prologue: Compute block 0
    if (num_blocks > 0) {
        compute_block(0, 0);
        MPI_Ireduce_scatter(send_buffers[0].data(), recv_buffers[0].data(), recvcounts.data(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD, &reqs[0]);
    }

    // Main loop
    for (int ib = 1; ib < num_blocks; ++ib) {
        int curr_buf = ib % 2;
        int prev_buf = (ib - 1) % 2;

        // Compute current block
        compute_block(ib, curr_buf);

        // Wait for previous comm to finish
        MPI_Wait(&reqs[prev_buf], MPI_STATUS_IGNORE);
        
        // Unpack previous block
        unpack_block(ib - 1, prev_buf);

        // Start current comm
        MPI_Ireduce_scatter(send_buffers[curr_buf].data(), recv_buffers[curr_buf].data(), recvcounts.data(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD, &reqs[curr_buf]);
    }

    // Epilogue: Wait and unpack last block
    if (num_blocks > 0) {
        int last_buf = (num_blocks - 1) % 2;
        MPI_Wait(&reqs[last_buf], MPI_STATUS_IGNORE);
        unpack_block(num_blocks - 1, last_buf);
    }

    print_distributed_matrix(A, local_rows, M, rank, size, "Final Matrix A");

    MPI_Finalize();
    return 0;
}
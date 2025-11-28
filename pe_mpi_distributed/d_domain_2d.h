#pragma once

#include <mpi.h>
#include <vector>
#include <iostream>

/**
 * @brief Distributed Domain 2D
 * 
 * Manages the decomposition of a 2D global domain into local slabs for MPI processing.
 * Currently supports 1D slab decomposition along the X-axis (rows).
 */
class DDomain2D
{
public:
    /**
     * @brief Construct a new DDomain2D object
     * 
     * @param in_comm MPI Communicator
     * @param in_global_nx Global size in X
     * @param in_global_ny Global size in Y
     * @param in_hx Grid spacing in X
     * @param in_hy Grid spacing in Y
     */
    DDomain2D(MPI_Comm in_comm, int in_global_nx, int in_global_ny, double in_hx, double in_hy);
    ~DDomain2D() = default;

    // Getters
    int get_global_nx() const { return global_nx; }
    int get_global_ny() const { return global_ny; }
    double get_hx() const { return hx; }
    double get_hy() const { return hy; }

    int get_rank() const { return rank; }
    int get_size() const { return size; }
    MPI_Comm get_comm() const { return comm; }

    // Local decomposition info (i-slab)
    int get_local_nx() const { return local_nx; }
    int get_local_ny() const { return local_ny; } // Usually same as global_ny for i-slab
    int get_local_i_start() const { return i_displs[rank]; }
    
    // Global decomposition info
    const std::vector<int>& get_i_counts() const { return i_counts; }
    const std::vector<int>& get_i_displs() const { return i_displs; }

    // Helper for J-slab decomposition (used for transpose)
    void get_j_decomposition(std::vector<int>& out_j_counts, std::vector<int>& out_j_displs) const;

private:
    MPI_Comm comm;
    int rank = 0;
    int size = 1;

    int global_nx = 0;
    int global_ny = 0;
    double hx = 0.0;
    double hy = 0.0;

    // I-slab decomposition (split X)
    std::vector<int> i_counts;
    std::vector<int> i_displs;
    int local_nx = 0;
    int local_ny = 0;

    static void split_1d(int n, int p, std::vector<int>& counts, std::vector<int>& displs);
};

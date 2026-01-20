#pragma once

#include "base/pch.h"

#include "base/domain/domain2d.h"
#include "base/domain/variable.h"
#include "base/location_boundary.h"
#include "io/config.h"
#include "pe/concat/domain_solver.h"
#include "pe/poisson/chasing_method2d.h"
#include "pe/poisson/poisson_fft2d.h"


#include <mpi.h>
#include <vector>

class MPIPoissonSolver2D : public DomainSolver2D
{
public:
    MPIPoissonSolver2D() {}

    MPIPoissonSolver2D(int                in_nx,
                       int                in_ny,
                       double             in_hx,
                       double             in_hy,
                       PDEBoundaryType    in_boundary_type_left,
                       PDEBoundaryType    in_boundary_type_right,
                       PDEBoundaryType    in_boundary_type_down,
                       PDEBoundaryType    in_boundary_type_up,
                       int                in_num_proc,
                       int                in_start_rank = 0,
                       EnvironmentConfig* in_env_config = nullptr,
                       MPI_Comm           in_comm       = MPI_COMM_WORLD);

    MPIPoissonSolver2D(Domain2DUniform*   in_domain,
                       Variable*          in_variable,
                       int                in_num_proc,
                       int                in_start_rank = 0,
                       EnvironmentConfig* in_env_config = nullptr,
                       MPI_Comm           in_comm       = MPI_COMM_WORLD);

    ~MPIPoissonSolver2D();

    void solve(field2& f, bool is_debugmode = true) override;
    void solve_collective_root_owned(field2& f, bool is_debugmode = true) override;
    bool is_comm_root() const override { return is_active && active_rank == 0; }

    double get_hx() const override { return hx; }
    double get_hy() const override { return hy; }

private:
    // Global problem size
    int    nx = 0, ny = 0;
    double hx = 0.0, hy = 0.0;

    // Domain/variable and env
    Variable*          var        = nullptr;
    Domain2DUniform*   domain     = nullptr;
    EnvironmentConfig* env_config = nullptr;

    // Boundary types
    PDEBoundaryType boundary_type_left  = PDEBoundaryType::Null;
    PDEBoundaryType boundary_type_right = PDEBoundaryType::Null;
    PDEBoundaryType boundary_type_down  = PDEBoundaryType::Null;
    PDEBoundaryType boundary_type_up    = PDEBoundaryType::Null;

    // MPI
    MPI_Comm world_comm  = MPI_COMM_WORLD;
    MPI_Comm active_comm = MPI_COMM_NULL;
    int      world_rank = 0, world_size = 1;
    int      active_rank = -1, active_size = 0;
    int      num_proc   = 1;
    int      start_rank = 0; // 起始世界进程编号（该进程成为子通信器的 rank 0）
    bool     is_active  = false;

    // Decomposition (i-slab for first stage, j-slab after transpose)
    std::vector<int> i_counts, i_displs; // split nx over active_size
    std::vector<int> j_counts, j_displs; // split ny over active_size

    // Local sizes
    int local_i_count = 0; // owned i count (x-direction slabs)
    int local_j_count = 0; // owned j count (y-direction slabs after transpose)

    // Local solvers
    PoissonFFT2D*    poisson_fft_y    = nullptr; // local y-direction FFT per i
    ChasingMethod2D* chasing_method_x = nullptr; // local x-direction chasing on transposed slabs

    // Local diagonal for x-chasing (subset for owned j indices)
    double* local_x_diag = nullptr;

    // Persistent work buffers (allocated after decomposition)
    field2 f_local;      // (local_i_count, ny)
    field2 fhat_local;   // (local_i_count, ny)
    field2 fhat_local_T; // (local_j_count, nx)
    field2 phat_local_T; // (local_j_count, nx)
    field2 phat_local;   // (local_i_count, ny)
    field2 p_local;      // (local_i_count, ny)

private:
    void setup_comm(MPI_Comm in_comm);
    void setup_problem_params();
    void build_decomposition();
    void build_local_solvers();
    void release_local_solvers();

    // Exact numeric helpers (must match OpenMP version)
    void cal_lambda(std::vector<double>& lambda_y_out) const;
    void build_local_x_diag();
    void boundary_assembly(field2& f);

    // Distributed ops
    void scatter_global_f_to_local(const field2& f_global, field2& f_local) const;
    void gather_local_to_global_f(const field2& f_local, field2& f_global) const;

    // Distributed transpose between i-slab (local_i_count x ny) and j-slab (local_j_count x nx)
    void distributed_transpose_i_to_j(const field2& in_i_slab, field2& out_j_slab) const;
    void distributed_transpose_j_to_i(const field2& in_j_slab, field2& out_i_slab) const;

    // Util
    static void split_1d(int n, int p, std::vector<int>& counts, std::vector<int>& displs);
};

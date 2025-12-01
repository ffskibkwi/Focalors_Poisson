#pragma once

#include "base/location_boundary.h"
#include "d_domain_2d.h"
#include "d_field_2d.h"
#include "io/config.h"
#include "pe/poisson/chasing_method2d.h"
#include "pe/poisson/poisson_fft2d.h"

#include <mpi.h>
#include <vector>

/**
 * @brief Distributed MPI Poisson Solver 2D
 *
 * Solves Poisson equation on a distributed domain (DDomain2D).
 * Input and Output are DField2D objects which hold local data.
 */
class MPIDistributedPoissonSolver2D
{
public:
    /**
     * @brief Construct a new MPIDistributedPoissonSolver2D object
     *
     * @param in_domain Pointer to the distributed domain structure
     * @param in_bc_left Boundary condition type Left
     * @param in_bc_right Boundary condition type Right
     * @param in_bc_down Boundary condition type Down
     * @param in_bc_up Boundary condition type Up
     * @param in_env_config Environment config (optional)
     */
    MPIDistributedPoissonSolver2D(DDomain2D*         in_domain,
                                  PDEBoundaryType    in_bc_left,
                                  PDEBoundaryType    in_bc_right,
                                  PDEBoundaryType    in_bc_down,
                                  PDEBoundaryType    in_bc_up,
                                  EnvironmentConfig* in_env_config = nullptr);

    ~MPIDistributedPoissonSolver2D();

    /**
     * @brief Solve the Poisson equation: Laplacian(u) = f
     *
     * @param f Input: RHS field (f). Output: Solution field (u).
     *          Data must be distributed in f before calling.
     */
    void solve(DField2D& f);

private:
    DDomain2D*         domain     = nullptr;
    EnvironmentConfig* env_config = nullptr;

    PDEBoundaryType bc_left;
    PDEBoundaryType bc_right;
    PDEBoundaryType bc_down;
    PDEBoundaryType bc_up;

    // Local solvers
    PoissonFFT2D*    poisson_fft_y    = nullptr;
    ChasingMethod2D* chasing_method_x = nullptr;

    // Local x-diagonal for chasing method
    double* local_x_diag = nullptr;

    // Work buffers
    field2 fhat_local;   // (local_i_count, ny)
    field2 fhat_local_T; // (local_j_count, nx)
    field2 phat_local_T; // (local_j_count, nx)
    field2 phat_local;   // (local_i_count, ny)

    // Helper methods
    void build_local_solvers();
    void release_local_solvers();
    void build_local_x_diag();
    void cal_lambda(std::vector<double>& lambda_y_out) const;

    // Boundary assembly on local parts
    void boundary_assembly_local(DField2D& f);

    // Distributed Transpose
    void distributed_transpose_i_to_j(const field2& in_i_slab, field2& out_j_slab) const;
    void distributed_transpose_j_to_i(const field2& in_j_slab, field2& out_i_slab) const;
};

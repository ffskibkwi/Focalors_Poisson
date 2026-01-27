#pragma once

#include "base/pch.h"

#include "base/domain/domain2d.h"
#include "base/domain/variable.h"
#include "io/config.h"
#include "poisson_solver2d_base.h"

#include <mpi.h>

class PoissonSolver2DSlabX : public PoissonSolver2DBase
{
public:
    PoissonSolver2DSlabX() {};
    PoissonSolver2DSlabX(Domain2DUniform*   in_domain,
                         Variable*          in_variable,
                         EnvironmentConfig* in_env_config   = nullptr,
                         MPI_Comm           in_communicator = MPI_COMM_WORLD);
    ~PoissonSolver2DSlabX();

    void solve(field2& f) override;

private:
    field2 f_hat, f_hat_T, p_hat_T, p_hat;

    std::string domain_name;

    EnvironmentConfig* env_config = nullptr;

    // slab
    MPI_Comm communicator;
    int      mpi_rank, mpi_size;
    int      snx, sny, ny_disp;
    double * buf_slab_x = nullptr, *buf_slab_y = nullptr;
};
#pragma once

#include "base/domain/variable2d_slab_x.h"
#include "concat_solver2d.h"

class ConcatPoissonSolver2DSlabX : public ConcatPoissonSolver2D
{
public:
    Variable2DSlabX* variable = nullptr;

    ConcatPoissonSolver2DSlabX(Variable2DSlabX* in_variable, EnvironmentConfig* in_env_config = nullptr)
        : ConcatPoissonSolver2D(in_variable, in_env_config)
        , variable(in_variable)
    {}

    void solve();

protected:
    std::unordered_map<Domain2DUniform*, DomainSolver2D*>
         construct_local_solver_map(Domain2DUniform* domain, bool is_local, MPI_Comm local_comm);
    void construct_solver_map_at_domain(Domain2DMPIUniform* domain);
    void construct_solver_map();

    void boundary_assembly();
};

#pragma once

#include "base/domain/variable2d_slab_x.h"
#include "concat_solver2d.h"

class ConcatPoissonSolver2DSlabX : public ConcatPoissonSolver2D
{
public:
    Variable2DSlabX* variable = nullptr;

    ConcatPoissonSolver2DSlabX(Variable2DSlabX* in_variable, EnvironmentConfig* in_env_config = nullptr);

    void solve();

protected:
    std::unordered_map<Domain2DUniform*, DomainSolver2D*>
         construct_local_solver_map(Domain2DUniform* domain, bool is_local, MPI_Comm local_comm);
    void construct_solver_map_at_domain(Domain2DUniformMPI* domain);
    void construct_solver_map();

    void boundary_assembly();

    void bond_add_slab(Domain2DUniformMPI*         domain_src,
                       Domain2DUniformMPI*         domain_dest,
                       LocationType                location,
                       double                      coeff,
                       field2*                     f_src,
                       const std::vector<field2*>& f_dest);

    double* get_buffer(int size);

    std::unordered_map<int, double*> buffer_map;
};

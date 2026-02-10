#pragma once

#include "base/pch.h"

#include "base/config.h"
#include "base/domain/domain2d.h"
#include "base/domain/variable2d.h"
#include "poisson_solver2d_base.h"

class PoissonSolver2D : public PoissonSolver2DBase
{
public:
    PoissonSolver2D() {};
    PoissonSolver2D(int             in_nx,
                    int             in_ny,
                    double          in_hx,
                    double          in_hy,
                    PDEBoundaryType in_boundary_type_left,
                    PDEBoundaryType in_boundary_type_right,
                    PDEBoundaryType in_boundary_type_down,
                    PDEBoundaryType in_boundary_type_up);
    PoissonSolver2D(Domain2DUniform* in_domain, Variable2D* in_variable = nullptr);

    void solve(field2& f) override;

private:
    void init();

    field2 buffer;

    std::string domain_name;

    EnvironmentConfig* env_config = nullptr;
};
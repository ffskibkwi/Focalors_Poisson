#pragma once

#include "base/pch.h"

#include "base/config.h"
#include "base/domain/domain3d.h"
#include "base/domain/variable3d.h"
#include "poisson_solver3d_base.h"

class PoissonSolver3D : public PoissonSolver3DBase
{
public:
    PoissonSolver3D() {};
    PoissonSolver3D(int             in_nx,
                    int             in_ny,
                    int             in_nz,
                    double          in_hx,
                    double          in_hy,
                    double          in_hz,
                    PDEBoundaryType in_boundary_type_xneg,
                    PDEBoundaryType in_boundary_type_xpos,
                    PDEBoundaryType in_boundary_type_yneg,
                    PDEBoundaryType in_boundary_type_ypos,
                    PDEBoundaryType in_boundary_type_zneg,
                    PDEBoundaryType in_boundary_type_zpos);
    PoissonSolver3D(Domain3DUniform* in_domain, Variable3D* in_variable = nullptr);

    void solve(field3& f) override;

private:
    void init();

    field3 buffer;

    std::string domain_name;

    EnvironmentConfig* env_config = nullptr;
};
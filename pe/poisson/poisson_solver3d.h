#pragma once

#include "base/pch.h"

#include "base/domain/domain3d.h"
#include "base/domain/variable3d.h"
#include "io/config.h"
#include "poisson_solver3d_base.h"

class PoissonSolver3D : public PoissonSolver3DBase
{
public:
    PoissonSolver3D() {};
    PoissonSolver3D(Domain3DUniform* in_domain, Variable3D* in_variable, EnvironmentConfig* in_env_config = nullptr);

    void solve(field3& f) override;

private:
    field3 buffer;

    std::string domain_name;

    EnvironmentConfig* env_config = nullptr;
};
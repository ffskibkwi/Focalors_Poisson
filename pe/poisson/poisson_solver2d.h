#pragma once

#include "base/pch.h"

#include "base/domain/domain2d.h"
#include "base/domain/variable.h"
#include "io/config.h"
#include "poisson_solver2d_base.h"

class PoissonSolver2D : public PoissonSolver2DBase
{
public:
    PoissonSolver2D() {};
    PoissonSolver2D(Domain2DUniform* in_domain, Variable* in_variable, EnvironmentConfig* in_env_config = nullptr);

    void solve(field2& f) override;

private:
    field2 buffer;

    std::string domain_name;

    EnvironmentConfig* env_config = nullptr;
};
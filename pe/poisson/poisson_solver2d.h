#pragma once

#include "base/pch.h"

#include "base/location_boundary.h"
#include "base/domain/domain2d.h"
#include "poisson_fft2d.h"
#include "chasing_method2d.h"
#include "pe/concat/domain_solver.h"
#include "io/config.h"
#include "base/domain/variable.h"

class PoissonSolver2D : public DomainSolver2D
{
public:
    PoissonSolver2D() {};

    PoissonSolver2D(int in_nx, int in_ny, double in_hx, double in_hy, PDEBoundaryType in_BoundaryTypeXNegative, PDEBoundaryType in_BoundaryTypeXPositive, PDEBoundaryType in_BoundaryTypeYNegative, PDEBoundaryType in_BoundaryTypeYPositive, EnvironmentConfig* in_env_config = nullptr);
    PoissonSolver2D(Domain2DUniform* in_domain, Variable* in_variable, EnvironmentConfig* in_env_config = nullptr);
    ~PoissonSolver2D();

    void init();

    void cal_lambda();

    void solve(field2& f) override;

private:
    EnvironmentConfig* env_config = nullptr;

private:
    int    nx, ny;
    double hx, hy;
    field2 buffer;

    PDEBoundaryType BoundaryTypeXNegative;
    PDEBoundaryType BoundaryTypeXPositive;
    PDEBoundaryType BoundaryTypeYNegative;
    PDEBoundaryType BoundaryTypeYPositive;

    double* lambda_y;
    double* x_diag; // The diagonal elements of the x direction subequation

    PoissonFFT2D*    poisson_fft_y;
    ChasingMethod2D* chasing_method_x;
};
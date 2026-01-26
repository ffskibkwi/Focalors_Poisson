#pragma once

#include "base/pch.h"

#include "base/domain/domain2d.h"
#include "base/domain/variable.h"
#include "base/location_boundary.h"
#include "chasing_method2d.h"
#include "io/config.h"
#include "pe/concat/domain_solver.h"
#include "poisson_fft2d.h"

class PoissonSolver2D : public DomainSolver2D
{
public:
    PoissonSolver2D() {};

    PoissonSolver2D(int                in_nx,
                    int                in_ny,
                    double             in_hx,
                    double             in_hy,
                    PDEBoundaryType    in_boundary_type_left,
                    PDEBoundaryType    in_boundary_type_right,
                    PDEBoundaryType    in_boundary_type_down,
                    PDEBoundaryType    in_boundary_type_up,
                    EnvironmentConfig* in_env_config = nullptr);
    PoissonSolver2D(Domain2DUniform* in_domain, Variable* in_variable, EnvironmentConfig* in_env_config = nullptr);
    ~PoissonSolver2D();

    void init();

    void solve(field2& f) override;

    double get_hx() const override { return hx; }
    double get_hy() const override { return hy; }

private:
    int    nx, ny;
    double hx, hy;
    field2 buffer;

    Variable*        var    = nullptr;
    Domain2DUniform* domain = nullptr;

    EnvironmentConfig* env_config = nullptr;

    PDEBoundaryType boundary_type_left;
    PDEBoundaryType boundary_type_right;
    PDEBoundaryType boundary_type_down;
    PDEBoundaryType boundary_type_up;

    double* lambda_y;
    double* x_diag; // The diagonal elements of the x direction subequation

    PoissonFFT2D*    poisson_fft_y;
    ChasingMethod2D* chasing_method_x;

    void cal_lambda();

    int solve_call_count = 0;
};
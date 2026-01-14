#pragma once

#include "base/pch.h"

#include "base/domain/domain3d.h"
#include "base/domain/variable.h"
#include "base/domain/variable3d.h"
#include "base/location_boundary.h"
#include "chasing_method3d.h"
#include "io/config.h"
#include "pe/concat/domain_solver.h"
#include "poisson_fft3d.h"

class PoissonSolver3D : public DomainSolver3D
{
public:
    PoissonSolver3D() {};

    PoissonSolver3D(int                in_nx,
                    int                in_ny,
                    int                in_nz,
                    double             in_hx,
                    double             in_hy,
                    double             in_hz,
                    PDEBoundaryType    in_boundary_type_left,
                    PDEBoundaryType    in_boundary_type_right,
                    PDEBoundaryType    in_boundary_type_front,
                    PDEBoundaryType    in_boundary_type_back,
                    PDEBoundaryType    in_boundary_type_down,
                    PDEBoundaryType    in_boundary_type_up,
                    EnvironmentConfig* in_env_config = nullptr);
    PoissonSolver3D(Domain3DUniform* in_domain, Variable3D* in_variable, EnvironmentConfig* in_env_config = nullptr);
    ~PoissonSolver3D();

    void init();

    void solve(field3& f) override;

private:
    int    nx, ny, nz;
    double hx, hy, hz;
    field3 buffer;

    Variable3D*        var    = nullptr;
    Domain3DUniform* domain = nullptr;

    EnvironmentConfig* env_config = nullptr;

    PDEBoundaryType boundary_type_left;
    PDEBoundaryType boundary_type_right;
    PDEBoundaryType boundary_type_front;
    PDEBoundaryType boundary_type_back;
    PDEBoundaryType boundary_type_down;
    PDEBoundaryType boundary_type_up;

    double** x_diag; // The diagonal elements of the x direction subequation

    PoissonFFT3D*    poisson_fft_z;
    PoissonFFT3D*    poisson_fft_y;
    ChasingMethod3D* chasing_method_x;

    void cal_lambda(double* lambda, int global_length, int begin, int local_length, PDEBoundaryType BoundTypeNegative, PDEBoundaryType BoundTypePositive);
    void boundary_assembly(field3& f);
};
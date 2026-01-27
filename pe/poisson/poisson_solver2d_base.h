#pragma once

#include "base/pch.h"

#include "base/location_boundary.h"
#include "chasing_method2d.h"
#include "pe/concat/domain_solver.h"
#include "poisson_fft2d.h"

class PoissonSolver2DBase : public DomainSolver2D
{
public:
    PoissonSolver2DBase() {};
    PoissonSolver2DBase(int             in_nx,
                        int             in_ny,
                        double          in_hx,
                        double          in_hy,
                        PDEBoundaryType in_boundary_type_left,
                        PDEBoundaryType in_boundary_type_right,
                        PDEBoundaryType in_boundary_type_down,
                        PDEBoundaryType in_boundary_type_up);
    ~PoissonSolver2DBase();

    virtual void solve(field2& f) {}

protected:
    int    nx, ny;
    double hx, hy;

    PDEBoundaryType boundary_type_left;
    PDEBoundaryType boundary_type_right;
    PDEBoundaryType boundary_type_down;
    PDEBoundaryType boundary_type_up;

    double* x_diag; // The diagonal elements of the x direction subequation

    PoissonFFT2D*    poisson_fft_y;
    ChasingMethod2D* chasing_method_x;

    // This function is calcualting the lambda of the equation with all the diagonal elements equal to 0
    void cal_lambda(double*         lambda,
                    int             global_length,
                    int             begin,
                    int             local_length,
                    PDEBoundaryType BoundTypeNegative,
                    PDEBoundaryType BoundTypePositive);

    static void create_fft(PoissonFFT2D*& fft, PDEBoundaryType bound_neg, PDEBoundaryType bound_pos, int n1, int n2);

    int solve_call_count = 0;
};
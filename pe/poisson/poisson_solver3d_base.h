#pragma once

#include "base/pch.h"

#include "base/location_boundary.h"
#include "chasing_method3d.h"
#include "pe/concat/domain_solver.h"
#include "poisson_fft3d.h"

class PoissonSolver3DBase : public DomainSolver3D
{
public:
    PoissonSolver3DBase() {};
    PoissonSolver3DBase(int             in_nx,
                        int             in_ny,
                        int             in_nz,
                        double          in_hx,
                        double          in_hy,
                        double          in_hz,
                        PDEBoundaryType in_boundary_type_left,
                        PDEBoundaryType in_boundary_type_right,
                        PDEBoundaryType in_boundary_type_front,
                        PDEBoundaryType in_boundary_type_back,
                        PDEBoundaryType in_boundary_type_down,
                        PDEBoundaryType in_boundary_type_up);
    ~PoissonSolver3DBase();

    virtual void solve(field3& f) {}

protected:
    int    nx, ny, nz;
    double hx, hy, hz;

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

    // This function is calcualting the lambda of the equation with all the diagonal elements equal to 0
    void cal_lambda(double*         lambda,
                    int             global_length,
                    int             begin,
                    int             local_length,
                    PDEBoundaryType BoundTypeNegative,
                    PDEBoundaryType BoundTypePositive);

    static void
    create_fft(PoissonFFT3D*& fft, PDEBoundaryType bound_neg, PDEBoundaryType bound_pos, int n1, int n2, int n3);
};
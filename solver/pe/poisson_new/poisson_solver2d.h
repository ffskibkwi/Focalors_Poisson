#pragma once

#include "pch.h"

#include "core/boundary/boundary_type.h"
#include "core/domain/domain2d.h"
#include "poisson_fft2d.h"
#include "chasing_method2d.h"
#include "domain_solver.h"

class PoissonSolver2D : public DomainSolver2D
{
public:
    PoissonSolver2D() {};

    PoissonSolver2D(int in_nx, int in_ny, double in_hx, double in_hy, PDEBoundaryType in_BoundaryTypeXNegative, PDEBoundaryType in_BoundaryTypeXPositive, PDEBoundaryType in_BoundaryTypeYNegative, PDEBoundaryType in_BoundaryTypeYPositive);
    PoissonSolver2D(Domain2DUniform& in_domain);
    ~PoissonSolver2D();

    void init() override;

    void cal_lambda();

    void solve(field2& f) override;

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

    /*
    // for MPI
    void cal_lambda(double*           lambda,
                    int               global_length,
                    int               begin,
                    int               local_length,
                    PDEBoundaryType BoundaryTypeNegative,
                    PDEBoundaryType BoundaryTypePositive)
    {
        // This function is calcualting the lambda of the equation with all the diagonal elements equal to 0
        for (int i = begin; i <= begin + local_length - 1; i++)
        {
            if (BoundaryTypeNegative == PDEBoundaryType::Periodic &&
                BoundaryTypePositive == PDEBoundaryType::Periodic) // P-P
            {
                lambda[i - begin] = -2.0 * std::cos(2.0 * pi / global_length * std::floor(i / 2.0));
            }
            else if (BoundaryTypeNegative == PDEBoundaryType::Neumann &&
                     BoundaryTypePositive == PDEBoundaryType::Neumann) // N-N
            {
                lambda[i - begin] = -2.0 * std::cos(pi / global_length * i);
            }
            else if (BoundaryTypeNegative == PDEBoundaryType::Dirichlet &&
                     BoundaryTypePositive == PDEBoundaryType::Dirichlet) // D-D
            {
                lambda[i - begin] = -2.0 * std::cos(pi / (global_length + 1.0) * i);
            }
            else if ((BoundaryTypeNegative == PDEBoundaryType::Dirichlet &&
                      BoundaryTypePositive == PDEBoundaryType::Neumann) ||
                     (BoundaryTypeNegative == PDEBoundaryType::Neumann &&
                      BoundaryTypePositive == PDEBoundaryType::Dirichlet)) // D-N or N-D
            {
                lambda[i - begin] = -2.0 * std::cos(2.0 * pi * i / (2.0 * global_length + 1.0));
            }
        }
    }
    */
};
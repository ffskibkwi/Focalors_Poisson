#pragma once

#include "pch.h"

#include "chasing_method_2d.hpp"
#include "core/base/location_boundary.h"
#include "core/mesh_profile/mesh_profile_2d_uniform.h"
#include "poisson_fft_2d.hpp"
#include "poisson_solver_interface.h"

template<PDEBoundaryType BoundTypeXNegative,
         PDEBoundaryType BoundTypeXPositive,
         PDEBoundaryType BoundTypeYNegative,
         PDEBoundaryType BoundTypeYPositive>
class PoissonSolver2D : public PoissonSolver2DInterface
{
public:
    PoissonSolver2D() {}

    PoissonSolver2D(int in_nx, int in_ny, double in_hx = 1.0, double in_hy = 1.0)
    {
        nx = in_nx;
        ny = in_ny;

        hx = in_hx;
        hy = in_hy;

        set_output_size(nx, ny);
        buffer.init(output_nx, output_ny, "buffer");

        poisson_fft_y.init(nx, ny);

        lambda_y = new double[ny];
        x_diag   = new double[ny];

        cal_lambda(lambda_y, ny, 1, ny, BoundTypeYNegative, BoundTypeYPositive);

        for (int j = 0; j < ny; j++)
        {
            x_diag[j] = -2.0 + hx * hx / hy / hy * (lambda_y[j] - 2.0);
        }
        delete[] lambda_y;

        bool is_no_Dirichlet =
            (BoundTypeXNegative != PDEBoundaryType::Dirichlet && BoundTypeXPositive != PDEBoundaryType::Dirichlet &&
             BoundTypeYNegative != PDEBoundaryType::Dirichlet && BoundTypeYPositive != PDEBoundaryType::Dirichlet);
        bool has_last_vector = true;
        chasing_method_x.init(nx, ny, x_diag, is_no_Dirichlet, has_last_vector);
    }

    PoissonSolver2D(MeshProfile2DUniform& in_mesh_profile)
    {
        nx = in_mesh_profile.nx;
        ny = in_mesh_profile.ny;

        hx = in_mesh_profile.hx;
        hy = in_mesh_profile.hy;

        set_output_size(nx, ny);
        buffer.init(output_nx, output_ny, "buffer");

        poisson_fft_y.init(nx, ny);

        // Calculate the lambda
        lambda_y = new double[ny];
        x_diag   = new double[ny];

        cal_lambda(lambda_y, ny, 1, ny, BoundTypeYNegative, BoundTypeYPositive);

        for (int j = 0; j < ny; j++)
        {
            x_diag[j] = -2.0 + hx * hx / hy / hy * (lambda_y[j] - 2.0);
        }
        delete[] lambda_y;

        bool is_no_Dirichlet =
            (BoundTypeXNegative != PDEBoundaryType::Dirichlet && BoundTypeXPositive != PDEBoundaryType::Dirichlet &&
             BoundTypeYNegative != PDEBoundaryType::Dirichlet && BoundTypeYPositive != PDEBoundaryType::Dirichlet);
        bool has_last_vector = true;
        chasing_method_x.init(nx, ny, x_diag, is_no_Dirichlet, has_last_vector);
    }

    ~PoissonSolver2D()
    {
        delete[] x_diag;
    }

    void solve(field2& f) override
    {
        buffer.set_size(nx, ny);
        poisson_fft_y.transform(f, buffer);

        // 逐列/逐行追赶求解
        chasing_method_x.chasing(buffer, f);

        buffer.set_size(nx, ny);
        poisson_fft_y.transform_transpose(f, buffer);

        std::swap(f, buffer);
    }

private:
    int    nx, ny;
    double hx, hy;
    field2 buffer;

    double * lambda_y;
    double* x_diag; // The diagonal elements of the x direction subequation

    PoissonFFT2D<BoundTypeYNegative, BoundTypeYPositive>    poisson_fft_y;
    ChasingMethod2D<BoundTypeXNegative, BoundTypeXPositive> chasing_method_x;

    void cal_lambda(double*           lambda,
                    int               global_length,
                    int               begin,
                    int               local_length,
                    PDEBoundaryType BoundTypeNegative,
                    PDEBoundaryType BoundTypePositive)
    {
        // This function is calcualting the lambda of the equation with all the diagonal elements equal to 0
        for (int i = begin; i <= begin + local_length - 1; i++)
        {
            if (BoundTypeNegative == PDEBoundaryType::Periodic &&
                BoundTypePositive == PDEBoundaryType::Periodic) // P-P
            {
                lambda[i - begin] = -2.0 * std::cos(2.0 * pi / global_length * std::floor(i / 2.0));
            }
            else if (BoundTypeNegative == PDEBoundaryType::Neumann &&
                     BoundTypePositive == PDEBoundaryType::Neumann) // N-N
            {
                lambda[i - begin] = -2.0 * std::cos(pi / global_length * i);
            }
            else if (BoundTypeNegative == PDEBoundaryType::Dirichlet &&
                     BoundTypePositive == PDEBoundaryType::Dirichlet) // D-D
            {
                lambda[i - begin] = -2.0 * std::cos(pi / (global_length + 1.0) * i);
            }
            else if ((BoundTypeNegative == PDEBoundaryType::Dirichlet &&
                      BoundTypePositive == PDEBoundaryType::Neumann) ||
                     (BoundTypeNegative == PDEBoundaryType::Neumann &&
                      BoundTypePositive == PDEBoundaryType::Dirichlet)) // D-N or N-D
            {
                lambda[i - begin] = -2.0 * std::cos(2.0 * pi * i / (2.0 * global_length + 1.0));
            }
        }
    }
};
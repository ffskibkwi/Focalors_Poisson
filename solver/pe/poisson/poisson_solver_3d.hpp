#pragma once

#include "pch.h"

#include "chasing_method_3d.hpp"
#include "core/boundary/boundary_type.h"
#include "core/mesh_profile/mesh_profile_3d_uniform.h"
#include "poisson_fft_3d.hpp"
#include "poisson_solver_interface.h"

template<FluidBoundaryType BoundTypeXNegative,
         FluidBoundaryType BoundTypeXPositive,
         FluidBoundaryType BoundTypeYNegative,
         FluidBoundaryType BoundTypeYPositive,
         FluidBoundaryType BoundTypeZNegative,
         FluidBoundaryType BoundTypeZPositive>
class PoissonSolver3D : public PoissonSolver3DInterface
{
public:
    PoissonSolver3D() {}

    PoissonSolver3D(MeshProfile3DUniform& in_mesh_profile)
    {
        nx = in_mesh_profile.nx;
        ny = in_mesh_profile.ny;
        nz = in_mesh_profile.nz;

        hx = in_mesh_profile.hx;
        hy = in_mesh_profile.hy;
        hz = in_mesh_profile.hz;

        set_output_size(nx, ny, nz);
        buffer.init(output_nx, output_ny, output_nz, "buffer");

        poisson_fft_z.init(nx, ny, nz);
        poisson_fft_y.init(nz, nx, ny);

        // Calculate the lambda
        lambda_y = new double[ny];
        lambda_z = new double[nz];
        x_diag   = new double*[nz];

        cal_lambda(lambda_y, ny, 1, ny, BoundTypeYNegative, BoundTypeYPositive);
        cal_lambda(lambda_z, nz, 1, nz, BoundTypeZNegative, BoundTypeZPositive);

        for (int k = 0; k < nz; k++)
        {
            x_diag[k] = new double[ny];
            for (int j = 0; j < ny; j++)
            {
                x_diag[k][j] = -2.0 + hx * hx / hy / hy * (lambda_y[j] - 2.0) + hx * hx / hz / hz * (lambda_z[k] - 2.0);
            }
        }
        delete[] lambda_y;
        delete[] lambda_z;

        bool is_no_Dirichlet =
            (BoundTypeXNegative != FluidBoundaryType::Dirichlet && BoundTypeXPositive != FluidBoundaryType::Dirichlet &&
             BoundTypeYNegative != FluidBoundaryType::Dirichlet && BoundTypeYPositive != FluidBoundaryType::Dirichlet &&
             BoundTypeZNegative != FluidBoundaryType::Dirichlet && BoundTypeZPositive != FluidBoundaryType::Dirichlet);
        bool has_last_vector = true;
        chasing_method_x.init(nz, ny, nx, x_diag, is_no_Dirichlet, has_last_vector);
    }

    ~PoissonSolver3D()
    {
        for (int k = 0; k < nz; k++)
        {
            delete[] x_diag[k];
        }
        delete[] x_diag;
    }

    void solve(field3& f) override
    {
        buffer.set_size(nx, ny, nz);
        poisson_fft_z.transform(f, buffer);

        f.set_size(nz, nx, ny);
        buffer.transpose(f, {2, 0, 1});

        buffer.set_size(nz, nx, ny);
        poisson_fft_y.transform(f, buffer);

        f.set_size(nz, ny, nx);
        buffer.transpose(f, {0, 2, 1});

        buffer.set_size(nz, ny, nx);
        chasing_method_x.chasing(f, buffer);

        f.set_size(nz, nx, ny);
        buffer.transpose(f, {0, 2, 1});

        buffer.set_size(nz, nx, ny);
        poisson_fft_y.transform_transpose(f, buffer);

        f.set_size(nx, ny, nz);
        buffer.transpose(f, {1, 2, 0});

        buffer.set_size(output_nx, output_ny, output_nz);
        poisson_fft_z.transform_transpose(f, buffer);

        std::swap(f, buffer);
    }

private:
    int    nx, ny, nz;
    double hx, hy, hz;
    field3 buffer;

    double * lambda_y, *lambda_z;
    double** x_diag; // The diagonal elements of the x direction subequation

    PoissonFFT3D<BoundTypeZNegative, BoundTypeZPositive>    poisson_fft_z;
    PoissonFFT3D<BoundTypeYNegative, BoundTypeYPositive>    poisson_fft_y;
    ChasingMethod3D<BoundTypeXNegative, BoundTypeXPositive> chasing_method_x;

    void cal_lambda(double*           lambda,
                    int               global_length,
                    int               begin,
                    int               local_length,
                    FluidBoundaryType BoundTypeNegative,
                    FluidBoundaryType BoundTypePositive)
    {
        // This function is calcualting the lambda of the equation with all the diagonal elements equal to 0
        for (int i = begin; i <= begin + local_length - 1; i++)
        {
            if (BoundTypeNegative == FluidBoundaryType::Periodic &&
                BoundTypePositive == FluidBoundaryType::Periodic) // P-P
            {
                lambda[i - begin] = -2.0 * std::cos(2.0 * pi / global_length * std::floor(i / 2.0));
            }
            else if (BoundTypeNegative == FluidBoundaryType::Neumann &&
                     BoundTypePositive == FluidBoundaryType::Neumann) // N-N
            {
                lambda[i - begin] = -2.0 * std::cos(pi / global_length * i);
            }
            else if (BoundTypeNegative == FluidBoundaryType::Dirichlet &&
                     BoundTypePositive == FluidBoundaryType::Dirichlet) // D-D
            {
                lambda[i - begin] = -2.0 * std::cos(pi / (global_length + 1.0) * i);
            }
            else if ((BoundTypeNegative == FluidBoundaryType::Dirichlet &&
                      BoundTypePositive == FluidBoundaryType::Neumann) ||
                     (BoundTypeNegative == FluidBoundaryType::Neumann &&
                      BoundTypePositive == FluidBoundaryType::Dirichlet)) // D-N or N-D
            {
                lambda[i - begin] = -2.0 * std::cos(2.0 * pi * i / (2.0 * global_length + 1.0));
            }
        }
    }
};
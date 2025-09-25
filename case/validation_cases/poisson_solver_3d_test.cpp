#include "core/boundary/boundary_type.h"
#include "pch.h"
#include "pe/poisson/poisson_solver_3d.hpp"
#include <cmath>
#include <iomanip>
#include <iostream>

/**
 * @brief Test the PoissonSolver3D class to solve the 3D Poisson equation
 * Boundary conditions: periodic in x-direction, Neumann in y and z directions
 */
int main(int argc, char* argv[])
{
    // Initialize grid parameters
    int    nx = 12;
    int    ny = 24;
    int    nz = 16;
    double Lx = 2.0 * pi; // Periodic boundary length is 2π
    double Ly = 1.0;
    double Lz = 2.0;

    std::cout << "Creating uniform grid..." << std::endl;

    // Create uniform grid
    MeshProfile3DUniform mesh_profile(nx, ny, nz, Lx, Ly, Lz);

    std::cout << "Creating Poisson solver..." << std::endl;

    // Create Poisson solver with specified boundary conditions
    // x-direction: periodic
    // y-direction: Neumann
    // z-direction: Neumann
    PoissonSolver3D<PDEBoundaryType::Periodic,
                    PDEBoundaryType::Periodic,
                    PDEBoundaryType::Periodic,
                    PDEBoundaryType::Periodic,
                    PDEBoundaryType::Neumann,
                    PDEBoundaryType::Neumann>
        poisson_solver(mesh_profile);

    std::cout << "Creating source field and solution fields..." << std::endl;

    // Create source field, solution field, and analytical solution field
    field3 source(nx, ny, nz, "source");
    field3 solution(nx, ny, nz, "solution");
    field3 analytical_solution(nx, ny, nz, "analytical_solution");

    // Grid spacing
    double hx = mesh_profile.hx;
    double hy = mesh_profile.hy;
    double hz = mesh_profile.hz;

    std::cout << "hx = " << hx << ", hy = " << hy << ", hz = " << hz << std::endl;

    std::cout << "Setting analytical solution and computing source term..." << std::endl;

    // Define the analytical solution
    // u(x,y,z) = sin(x) * cos(y) * cos(z)
    for (int i = 0; i < nx; i++)
    {
        for (int j = 0; j < ny; j++)
        {
            for (int k = 0; k < nz; k++)
            {
                analytical_solution(i, j, k) = exp(0.3 * pi * (j + 0.5) * hy / Ly) *
                                               exp(0.2 * pi * (i + 0.5) * hx / Lx) *
                                               exp(2.0 * pi * (k + 0.5) * hz / Lz) * cos(2 * pi * (i + 0.5) * hx / Lx) *
                                               cos(2 * pi * (j + 0.5) * hy / Ly) * cos(2 * pi * (k + 0.5) * hz / Lz);
            }
        }
    }

    // Print the analytical solution values
    for (int i = 0; i < nx; i++)
    {
        for (int j = 0; j < ny; j++)
        {
            for (int k = 0; k < nz; k++)
            {
                std::cout << "analytical_solution[" << i << "][" << j << "][" << k
                          << "] = " << analytical_solution(i, j, k) << std::endl;
            }
        }
    }

    // Compute source term using finite difference Laplacian
    for (int i = 0; i < nx; i++)
    {
        for (int j = 0; j < ny; j++)
        {
            for (int k = 0; k < nz; k++)
            {
                // Compute Laplacian: ∇²u = ∂²u/∂x² + ∂²u/∂y² + ∂²u/∂z²

                // ∂²u/∂x² (Periodic)
                double d2udx2;
                if (i == 0)
                {
                    d2udx2 = (analytical_solution(1, j, k) - 2 * analytical_solution(0, j, k) +
                              analytical_solution(nx - 1, j, k)) /
                             (hx * hx);
                }
                else if (i == nx - 1)
                {
                    d2udx2 = (analytical_solution(0, j, k) - 2 * analytical_solution(nx - 1, j, k) +
                              analytical_solution(nx - 2, j, k)) /
                             (hx * hx);
                }
                else
                {
                    d2udx2 = (analytical_solution(i + 1, j, k) - 2 * analytical_solution(i, j, k) +
                              analytical_solution(i - 1, j, k)) /
                             (hx * hx);
                }

                // ∂²u/∂y² (Periodic — intentionally inconsistent with Neumann setting)
                double d2udy2;
                if (j == 0)
                {
                    d2udy2 = (analytical_solution(i, 1, k) - 2 * analytical_solution(i, 0, k) +
                              analytical_solution(i, ny - 1, k)) /
                             (hy * hy);
                }
                else if (j == ny - 1)
                {
                    d2udy2 = (analytical_solution(i, ny - 2, k) - 2 * analytical_solution(i, ny - 1, k) +
                              analytical_solution(i, 0, k)) /
                             (hy * hy);
                }
                else
                {
                    d2udy2 = (analytical_solution(i, j + 1, k) - 2 * analytical_solution(i, j, k) +
                              analytical_solution(i, j - 1, k)) /
                             (hy * hy);
                }

                // ∂²u/∂z² (Neumann)
                double d2udz2;
                if (k == 0)
                {
                    d2udz2 = (analytical_solution(i, j, 1) - analytical_solution(i, j, 0)) / (hz * hz);
                }
                else if (k == nz - 1)
                {
                    d2udz2 = (analytical_solution(i, j, nz - 2) - analytical_solution(i, j, nz - 1)) / (hz * hz);
                }
                else
                {
                    d2udz2 = (analytical_solution(i, j, k + 1) - 2 * analytical_solution(i, j, k) +
                              analytical_solution(i, j, k - 1)) /
                             (hz * hz);
                }

                // Source term is negative Laplacian
                source(i, j, k) = hx * hx * (d2udx2 + d2udy2 + d2udz2);
            }
        }
    }

    // Subtract mean of the source term
    double source_mean = source.mean();
    std::cout << "Source mean: " << source_mean << std::endl;

    for (int i = 0; i < nx; i++)
    {
        for (int j = 0; j < ny; j++)
        {
            for (int k = 0; k < nz; k++)
            {
                source(i, j, k) -= source_mean;
            }
        }
    }

    // Subtract mean of the analytical solution
    double analytical_mean = analytical_solution.mean();
    std::cout << "Analytical solution mean: " << analytical_mean << std::endl;

    for (int i = 0; i < nx; i++)
    {
        for (int j = 0; j < ny; j++)
        {
            for (int k = 0; k < nz; k++)
            {
                analytical_solution(i, j, k) -= analytical_mean;
            }
        }
    }

    std::cout << "Solving Poisson equation..." << std::endl;

    // Solve Poisson equation
    solution = source; // Copy source to solution field; solver operates in-place
    poisson_solver.solve(solution);

    // Subtract mean from computed solution
    double mean = solution.mean();
    std::cout << "Solution mean: " << mean << std::endl;

    for (int i = 0; i < nx; i++)
    {
        for (int j = 0; j < ny; j++)
        {
            for (int k = 0; k < nz; k++)
            {
                solution(i, j, k) -= mean;
            }
        }
    }

    std::cout << "Solving Poisson equation..." << std::endl;

    // Compute error
    double max_error = 0.0;
    double l2_error  = 0.0;

    for (int i = 0; i < nx; i++)
    {
        for (int j = 0; j < ny; j++)
        {
            for (int k = 0; k < nz; k++)
            {
                double error = std::abs(solution(i, j, k) - analytical_solution(i, j, k));
                max_error    = std::max(max_error, error);
                l2_error += error * error;
                if (error * error > 1e-10)
                {
                    std::cout << "x = " << i << ", y = " << j << ", z = " << k << ", error = " << error << std::endl;
                }
            }
        }
    }

    l2_error = std::sqrt(l2_error / (nx * ny * nz));

    std::cout << "Solving Poisson equation..." << std::endl;

    // Print result on a yz-plane at a fixed x index
    int x_index = nx / 4;

    std::cout << "\nResults on yz-plane at x = " << x_index << " (x = " << x_index * hx << "):" << std::endl;

    std::cout << "\nNumerical solution matrix:" << std::endl;
    std::cout << std::setw(10) << "y\\z";
    for (int k = 0; k < nz; k++)
        std::cout << std::setw(10) << k;
    std::cout << std::endl;

    for (int j = 0; j < ny; j++)
    {
        std::cout << std::setw(10) << j;
        for (int k = 0; k < nz; k++)
        {
            std::cout << std::setw(10) << std::fixed << std::setprecision(4) << solution(x_index, j, k);
        }
        std::cout << std::endl;
    }

    std::cout << "\nAnalytical solution matrix:" << std::endl;
    std::cout << std::setw(10) << "y\\z";
    for (int k = 0; k < nz; k++)
        std::cout << std::setw(10) << k;
    std::cout << std::endl;

    for (int j = 0; j < ny; j++)
    {
        std::cout << std::setw(10) << j;
        for (int k = 0; k < nz; k++)
        {
            std::cout << std::setw(10) << std::fixed << std::setprecision(4) << analytical_solution(x_index, j, k);
        }
        std::cout << std::endl;
    }

    std::cout << "\nError matrix:" << std::endl;
    std::cout << std::setw(10) << "y\\z";
    for (int k = 0; k < nz; k++)
        std::cout << std::setw(10) << k;
    std::cout << std::endl;

    for (int j = 0; j < ny; j++)
    {
        std::cout << std::setw(10) << j;
        for (int k = 0; k < nz; k++)
        {
            std::cout << std::setw(10) << std::fixed << std::setprecision(4)
                      << std::abs(solution(x_index, j, k) - analytical_solution(x_index, j, k));
        }
        std::cout << std::endl;
    }

    std::cout << "\nFull comparison:" << std::endl;
    for (int i = 0; i < nx; i++)
    {
        std::cout << std::setw(10) << i;
        for (int j = 0; j < ny; j++)
        {
            for (int k = 0; k < nz; k++)
            {
                std::cout << "x = " << i << ", y = " << j << ", z = " << k << ", solution = " << solution(i, j, k)
                          << ", analytical_solution = " << analytical_solution(i, j, k)
                          << ", error = " << std::abs(solution(i, j, k) - analytical_solution(i, j, k)) << std::endl;
            }
        }
        std::cout << std::endl;
    }

    std::cout << "\nOverall error statistics:" << std::endl;
    std::cout << "Max error: " << max_error << std::endl;
    std::cout << "L2 error: " << l2_error << std::endl;

    return 0;
}

#include "pch.h"

#include "core/base/field2.h"
#include "core/base/location_boundary.h"
#include "pe/poisson/chasing_method_base.h"
#include <cmath>
#include <iomanip>
#include <iostream>
#include <mpi.h>
#include <vector>

/**
 * @brief Test program for the Thomas algorithm (Chasing method)
 *
 * This program tests various methods implemented in the ChasingMethodBase class for solving tridiagonal systems.
 * Tests include:
 * 1. Standard Thomas algorithm (without precomputation)
 * 2. Standard Thomas algorithm (with precomputation)
 * 3. Thomas algorithm with periodic boundary conditions
 * 4. Thomas algorithm for singular matrices
 */
class ChasingMethodTest : public ChasingMethodBase
{
public:
    // Test using protected methods
    void test_standard_chasing_without_precompute()
    {
        std::cout << "=== Testing Standard Thomas Algorithm (Without Precomputation) ===" << std::endl;

        const int n             = 10;
        double    x_diag_single = -2.1; // Diagonal element value

        std::vector<double> f(n, 0.0);     // Right-hand side vector
        std::vector<double> c(n - 1, 0.0); // Helper vector
        std::vector<double> y(n, 0.0);     // Helper vector
        std::vector<double> x(n, 0.0);     // Solution vector
        std::vector<double> exact(n, 0.0); // Exact solution

        // Set exact solution and compute RHS
        for (int i = 0; i < n; i++)
            exact[i] = sin(M_PI * i / (n - 1));
        for (int i = 0; i < n; i++)
        {
            if (i == 0)
                f[i] = x_diag_single * exact[i] + exact[i + 1];
            else if (i == n - 1)
                f[i] = x_diag_single * exact[i] + exact[i - 1];
            else
                f[i] = exact[i - 1] + x_diag_single * exact[i] + exact[i + 1];
        }

        std::cout << "Right-hand side vector f:" << std::endl;
        for (int i = 0; i < n; i++)
            std::cout << "f[" << i << "] = " << std::setprecision(6) << f[i] << std::endl;
        std::cout << std::endl;

        // Solve without precomputation
        chasing_standard_without_precompute(x_diag_single,
                                            n,
                                            f.data(),
                                            c.data(),
                                            y.data(),
                                            x.data(),
                                            PDEBoundaryType::Dirichlet,
                                            PDEBoundaryType::Dirichlet);

        // Compute and print error
        double error = 0.0;
        std::cout << "Index\tComputed\tExact\t\tError" << std::endl;
        for (int i = 0; i < n; i++)
        {
            double diff = std::abs(x[i] - exact[i]);
            error       = std::max(error, diff);
            std::cout << i << "\t" << std::setprecision(6) << x[i] << "\t\t" << exact[i] << "\t\t" << diff << std::endl;
        }
        std::cout << "Max Error: " << error << std::endl;
    }

    void test_standard_chasing_with_precompute()
    {
        std::cout << "\n=== Testing Standard Thomas Algorithm (With Precomputation) ===" << std::endl;

        const int n             = 10;
        double    x_diag_single = -2.1;

        std::vector<double> f(n, 0.0);     // Right-hand side vector
        std::vector<double> c(n - 1, 0.0); // Helper vector
        std::vector<double> y(n, 0.0);     // Helper vector
        std::vector<double> x(n, 0.0);     // Solution vector
        std::vector<double> exact(n, 0.0); // Exact solution

        // Precompute c vector
        chasing_standard_precompute(
            x_diag_single, n, c.data(), PDEBoundaryType::Dirichlet, PDEBoundaryType::Dirichlet);

        std::cout << "Helper vector c:" << std::endl;
        for (int i = 0; i < n - 1; i++)
            std::cout << "c[" << i << "] = " << std::setprecision(6) << c[i] << std::endl;
        std::cout << std::endl;

        for (int i = 0; i < n; i++)
            exact[i] = sin(M_PI * i / (n - 1));
        for (int i = 0; i < n; i++)
        {
            if (i == 0)
                f[i] = x_diag_single * exact[i] + exact[i + 1];
            else if (i == n - 1)
                f[i] = x_diag_single * exact[i] + exact[i - 1];
            else
                f[i] = exact[i - 1] + x_diag_single * exact[i] + exact[i + 1];
        }

        std::cout << "Right-hand side vector f:" << std::endl;
        for (int i = 0; i < n; i++)
            std::cout << "f[" << i << "] = " << std::setprecision(6) << f[i] << std::endl;
        std::cout << std::endl;

        // Solve with precomputation
        chasing_standard_with_precompute(x_diag_single,
                                         n,
                                         f.data(),
                                         c.data(),
                                         y.data(),
                                         x.data(),
                                         PDEBoundaryType::Dirichlet,
                                         PDEBoundaryType::Dirichlet);

        double error = 0.0;
        std::cout << "Index\tComputed\tExact\t\tError" << std::endl;
        for (int i = 0; i < n; i++)
        {
            double diff = std::abs(x[i] - exact[i]);
            error       = std::max(error, diff);
            std::cout << i << "\t" << std::setprecision(6) << x[i] << "\t\t" << exact[i] << "\t\t" << diff << std::endl;
        }
        std::cout << "Max Error: " << error << std::endl;
    }

    void test_periodic_chasing()
    {
        std::cout << "\n=== Testing Periodic Boundary Condition Thomas Algorithm ===" << std::endl;

        const int n             = 10;
        double    x_diag_single = -2.1;

        std::vector<double> f(n, 0.0);          // right hand side vector
        std::vector<double> c_diri(n - 1, 0.0); // Helper vector for Dirichlet boundary conditions
        std::vector<double> y_diri(n, 0.0);     // Helper vector for Dirichlet boundary conditions
        std::vector<double> z_1(n, 0.0);        // Helper vector for periodic boundary conditions
        std::vector<double> z_2(n, 0.0);        // Helper vector for periodic boundary conditions
        std::vector<double> y(n, 0.0);          // Helper vector
        std::vector<double> x(n, 0.0);          // Solution vector
        std::vector<double> exact(n, 0.0);      // Exact solution

        for (int i = 0; i < n; i++)
            exact[i] = sin(M_PI * i / (n - 1));
        for (int i = 0; i < n; i++)
        {
            if (i == 0)
                f[i] = exact[n - 1] + x_diag_single * exact[i] + exact[i + 1];
            else if (i == n - 1)
                f[i] = exact[i - 1] + x_diag_single * exact[i] + exact[0];
            else
                f[i] = exact[i - 1] + x_diag_single * exact[i] + exact[i + 1];
        }

        std::cout << "Right-hand side vector f:" << std::endl;
        for (int i = 0; i < n; i++)
            std::cout << "f[" << i << "] = " << std::setprecision(6) << f[i] << std::endl;
        std::cout << std::endl;

        // Precompute helper vectors for periodic case
        chasing_periodic_precompute(x_diag_single, n, z_1.data(), z_2.data(), c_diri.data());

        // Solve with periodic boundary conditions
        chasing_periodic_with_precompute(
            x_diag_single, n, f.data(), c_diri.data(), y_diri.data(), z_1.data(), z_2.data(), y.data(), x.data());

        double error = 0.0;
        std::cout << "Index\tComputed\tExact\t\tError" << std::endl;
        for (int i = 0; i < n; i++)
        {
            double diff = std::abs(x[i] - exact[i]);
            error       = std::max(error, diff);
            std::cout << i << "\t" << std::setprecision(6) << x[i] << "\t\t" << exact[i] << "\t\t" << diff << std::endl;
        }
        std::cout << "Max Error: " << error << std::endl;
    }

    void test_singular_chasing()
    {
        std::cout << "\n=== Testing Singular Matrix Thomas Algorithm ===" << std::endl;

        const int n = 10;

        std::vector<double> f(n, 0.0);
        std::vector<double> x(n, 0.0);

        // Ensure RHS satisfies solvability condition for singular matrix
        double sum = 0.0;
        for (int i = 0; i < n - 1; i++)
        {
            f[i] = i + 1.0;
            sum += f[i];
        }
        f[n - 1] = -sum;

        std::cout << "Right-hand side vector f:" << std::endl;
        for (int i = 0; i < n; i++)
            std::cout << "f[" << i << "] = " << std::setprecision(6) << f[i] << std::endl;
        std::cout << std::endl;

        chasing_standard_singular(n, f.data(), x.data());

        std::cout << "Index\tComputed" << std::endl;
        for (int i = 0; i < n; i++)
            std::cout << i << "\t" << std::setprecision(6) << x[i] << std::endl;
    }

    void test_periodic_singular_chasing()
    {
        std::cout << "\n=== Testing Periodic Singular Matrix Thomas Algorithm ===" << std::endl;

        const int n = 10;

        std::vector<double> f(n, 0.0);
        std::vector<double> x(n, 0.0);

        // Ensure RHS satisfies solvability condition for periodic singular matrix
        double sum = 0.0;
        for (int i = 0; i < n; i++)
        {
            f[i] = sin(2.0 * M_PI * i / n);
            sum += f[i];
        }
        if (std::abs(sum) > 1e-10)
        {
            for (int i = 0; i < n; i++)
                f[i] -= sum / n;
        }

        std::cout << "Right-hand side vector f:" << std::endl;
        for (int i = 0; i < n; i++)
            std::cout << "f[" << i << "] = " << std::setprecision(6) << f[i] << std::endl;
        std::cout << std::endl;

        chasing_periodic_singular(n, f.data(), x.data());

        std::cout << "Index\tComputed" << std::endl;
        for (int i = 0; i < n; i++)
            std::cout << i << "\t" << std::setprecision(6) << x[i] << std::endl;
    }
};

int main(int argc, char* argv[])
{
    MPI_Init(&argc, &argv);

    ChasingMethodTest test;

    test.test_standard_chasing_without_precompute();
    test.test_standard_chasing_with_precompute();
    test.test_periodic_chasing();
    test.test_singular_chasing();
    test.test_periodic_singular_chasing();

    MPI_Finalize();

    return 0;
}

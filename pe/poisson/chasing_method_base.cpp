#include "chasing_method_base.h"
#include "base/location_boundary.h"

/**
 * @brief The standard chasing method for tridiagonal equations without precomputation.
 *
 * This implementation solves a tridiagonal system of equations without precomputing
 * helper vectors.  It's suitable when the tridiagonal matrix changes frequently.
 *
 * @param x_diag_single The value of the main diagonal.  Assumed to be constant.
 * @param f            The right-hand side vector.
 * @param n            The size of the matrix (number of rows/columns).
 * @param c            The helper vector (size N - 1).  Must be pre-initialized. Modified during execution.
 * @param y            The helper vector (size N).  Modified during execution.
 * @param x            The solution vector.  Output.
 * @param BoundaryTypeStart The boundary type at the start (first row).
 * @param BoundaryTypeEnd   The boundary type at the end (last row).
 *
 * @note
 * The key difference between the precompute and non-precompute versions is that
 * the precompute version calculates the helper vector 'c' beforehand (typically
 * when the tridiagonal matrix remains constant across multiple solves). This version calculates 'c' inside.
 */
void ChasingMethodBase::chasing_standard_without_precompute(double            x_diag_single,
                                                            int               n,
                                                            double*           f,
                                                            double*           c,
                                                            double*           y,
                                                            double*           x,
                                                            PDEBoundaryType BoundaryTypeStart,
                                                            PDEBoundaryType BoundaryTypeEnd)
{
    c[0] = (BoundaryTypeStart == PDEBoundaryType::Neumann) ? 1.0 / (x_diag_single + 1.0) : 1.0 / x_diag_single;
    for (int i = 1; i < n - 1; i++)
    {
        c[i] = 1.0 / (x_diag_single - c[i - 1]);
    }

    y[0] = (BoundaryTypeStart == PDEBoundaryType::Neumann) ? f[0] / (x_diag_single + 1.0) : f[0] / x_diag_single;
    for (int i = 1; i < n - 1; i++)
    {
        y[i] = (f[i] - y[i - 1]) / (x_diag_single - c[i - 1]);
    }
    y[n - 1] = (BoundaryTypeEnd == PDEBoundaryType::Neumann) ? (f[n - 1] - y[n - 2]) / (x_diag_single + 1.0 - c[n - 2]) :
                                                              (f[n - 1] - y[n - 2]) / (x_diag_single - c[n - 2]);

    x[n - 1] = y[n - 1];
    for (int i = n - 2; i >= 0; i--)
    {
        x[i] = y[i] - c[i] * x[i + 1];
    }
}

/**
 * @brief The standard chasing method for tridiagonal equations with precomputed helper vector 'c'.
 *
 * This implementation solves a tridiagonal system of equations with the helper
 * vector 'c' already precomputed. This is advantageous when the tridiagonal
 * matrix remains constant across multiple solves, as it avoids redundant calculations.
 *
 * @param x_diag_single The value of the main diagonal. Assumed to be constant.
 * @param f            The right-hand side vector.
 * @param n            The size of the matrix (number of rows/columns).
 * @param c            The precomputed helper vector (size N-1).  Must be pre-calculated. Not modified during execution.
 * @param y            The helper vector (size N).  Modified during execution.
 * @param x            The solution vector. Output.
 * @param BoundaryTypeStart The boundary type at the start (first row).
 * @param BoundaryTypeEnd   The boundary type at the end (last row).
 *
 * @note
 *  In the precompute version, the helper vector 'c' is calculated and stored
 *  beforehand. This is typically done when the tridiagonal matrix remains constant
 *  across multiple solves. The non-precompute version calculates 'c' inside the
 *  solver, which is more suitable when the matrix changes frequently.
 */
void ChasingMethodBase::chasing_standard_with_precompute(double            x_diag_single,
                                                         int               n,
                                                         double*           f,
                                                         double*           c,
                                                         double*           y,
                                                         double*           x,
                                                         PDEBoundaryType BoundaryTypeStart,
                                                         PDEBoundaryType BoundaryTypeEnd)
{
    y[0] = (BoundaryTypeStart == PDEBoundaryType::Neumann) ? f[0] / (x_diag_single + 1.0) : f[0] / x_diag_single;
    for (int i = 1; i < n - 1; i++)
    {
        y[i] = (f[i] - y[i - 1]) / (x_diag_single - c[i - 1]);
    }
    y[n - 1] = (BoundaryTypeEnd == PDEBoundaryType::Neumann) ? (f[n - 1] - y[n - 2]) / (x_diag_single + 1.0 - c[n - 2]) :
                                                              (f[n - 1] - y[n - 2]) / (x_diag_single - c[n - 2]);

    x[n - 1] = y[n - 1];
    for (int i = n - 2; i >= 0; i--)
    {
        x[i] = y[i] - c[i] * x[i + 1];
    }
}

/**
 * @brief Precomputes the helper vector 'c' for the chasing_standard_with_precompute method.
 *
 * @param x_diag_single The value of the main diagonal. Assumed to be constant.
 * @param n            The size of the matrix (number of rows/columns).
 * @param c            The helper vector (size N-1) to be precomputed. Output.
 * @param BoundaryTypeStart The boundary type at the start (first row).
 * @param BoundaryTypeEnd   The boundary type at the end (last row).
 *
 * @note
 *  This function precalculates the 'c' vector, which is used as a helper in the
 *  `chasing_standard_with_precompute` method for solving tridiagonal systems of equations.
 *  This precomputation is beneficial when the tridiagonal matrix remains constant and
 *  is solved multiple times.
 */
void ChasingMethodBase::chasing_standard_precompute(double            x_diag_single,
                                                    int               n,
                                                    double*           c,
                                                    PDEBoundaryType BoundaryTypeStart,
                                                    PDEBoundaryType BoundaryTypeEnd)
{
    c[0] = (BoundaryTypeStart == PDEBoundaryType::Neumann) ? 1.0 / (x_diag_single + 1.0) : 1.0 / x_diag_single;
    for (int i = 1; i < n - 1; i++)
    {
        c[i] = 1.0 / (x_diag_single - c[i - 1]);
    }
}

/**
 * @brief The standard chasing method for tridiagonal equations with singular matrix.
 *
 * This implementation solves a tridiagonal system of equations with a singular matrix.
 * It handles the special case where the matrix is singular by using a constant value
 * and a recursive approach to compute the solution.
 *
 * @param n The size of the matrix (number of rows/columns).
 * @param f The right-hand side vector.
 * @param x The solution vector. Output.
 *
 * @note
 * The solution starts with a constant value p_const and builds the solution
 * recursively. This approach is specifically designed for singular matrices
 * where standard elimination methods would fail.
 */
void ChasingMethodBase::chasing_standard_singular(int n, double* f, double* x)
{
    x[0] = p_const;
    x[1] = p_const + f[0];
    for (int i = 2; i < n; i++)
    {
        x[i] = f[i - 1] - x[i - 2] + 2.0 * x[i - 1];
    }
}

/**
 * @brief The periodic chasing method for tridiagonal equations with precomputed helper vectors.
 *
 * This implementation solves a tridiagonal system of equations with periodic boundary conditions.
 * It uses precomputed helper vectors to efficiently solve the system by transforming the periodic
 * problem into a standard tridiagonal problem with corrections.
 *
 * @param x_diag_single The value of the main diagonal. Assumed to be constant.
 * @param n            The size of the matrix (number of rows/columns).
 * @param f            The right-hand side vector.
 * @param c_diri       The precomputed helper vector for Dirichlet boundary conditions.
 * @param y_diri       The helper vector for intermediate calculations.
 * @param z_1          The first precomputed special solution vector.
 * @param z_2          The second precomputed special solution vector.
 * @param y            The intermediate solution vector.
 * @param x            The final solution vector. Output.
 *
 * @note
 * This method uses the Sherman-Morrison formula to handle the periodic boundary conditions.
 * It first solves a standard tridiagonal system with Dirichlet boundary conditions, then
 * applies corrections using the precomputed special solutions z_1 and z_2.
 */
void ChasingMethodBase::chasing_periodic_with_precompute(double  x_diag_single,
                                                         int     n,
                                                         double* f,
                                                         double* c_diri,
                                                         double* y_diri,
                                                         double* z_1,
                                                         double* z_2,
                                                         double* y,
                                                         double* x)
{
    chasing_standard_with_precompute(
        x_diag_single, n, f, c_diri, y_diri, y, PDEBoundaryType::Dirichlet, PDEBoundaryType::Dirichlet);
    double gamma_2 =
        (y[n - 1] / (1.0 + z_1[n - 1]) - y[0] / z_1[0]) / (z_2[n - 1] / (1.0 + z_1[n - 1]) - (1.0 + z_2[0]) / z_1[0]);
    double gamma_1 = y[0] / z_1[0] - gamma_2 * ((1.0 + z_2[0]) / z_1[0]);
    for (int i = 0; i < n; i++)
    {
        x[i] = y[i] - gamma_1 * z_1[i] - gamma_2 * z_2[i];
    }
}

/**
 * @brief Precomputes the helper vectors for the chasing_periodic_with_precompute method.
 *
 * This function calculates the special solution vectors z_1 and z_2, and the helper vector c_diri,
 * which are used in the periodic chasing method. These vectors need to be computed only once
 * when the tridiagonal matrix remains constant across multiple solves.
 *
 * @param x_diag_single The value of the main diagonal. Assumed to be constant.
 * @param n            The size of the matrix (number of rows/columns).
 * @param z_1          The first special solution vector to be computed. Output.
 * @param z_2          The second special solution vector to be computed. Output.
 * @param c_diri       The helper vector for Dirichlet boundary conditions. Output.
 *
 * @note
 * This precomputation involves solving two standard tridiagonal systems with special
 * right-hand sides to obtain the z_1 and z_2 vectors, which represent the influence
 * of the periodic boundary conditions. The c_diri vector is computed for use in the
 * standard chasing method with Dirichlet boundary conditions.
 */
void ChasingMethodBase::chasing_periodic_precompute(double  x_diag_single,
                                                    int     n,
                                                    double* z_1,
                                                    double* z_2,
                                                    double* c_diri)
{
    // Usually precompute only run once, thus using the temporary pointer is acceptable
    double* c_temp  = new double[n - 1];
    double* y_temp  = new double[n];
    double* u1_temp = new double[n];
    double* u2_temp = new double[n];
    for (int i = 0; i < n; i++)
    {
        u1_temp[i] = (i == 0) ? 1.0 : 0.0;
        u2_temp[i] = (i == n - 1) ? 1.0 : 0.0;
    }
    chasing_standard_without_precompute(
        x_diag_single, n, u1_temp, c_temp, y_temp, z_1, PDEBoundaryType::Dirichlet, PDEBoundaryType::Dirichlet);
    chasing_standard_without_precompute(
        x_diag_single, n, u2_temp, c_temp, y_temp, z_2, PDEBoundaryType::Dirichlet, PDEBoundaryType::Dirichlet);

    delete[] c_temp;
    delete[] y_temp;
    delete[] u1_temp;
    delete[] u2_temp;

    c_diri[0] = 1.0 / x_diag_single;
    for (int i = 1; i < n - 1; i++)
    {
        c_diri[i] = 1.0 / (x_diag_single - c_diri[i - 1]);
    }
}

/**
 * @brief The periodic chasing method for tridiagonal equations with singular matrix.
 *
 * This implementation solves a tridiagonal system of equations with periodic boundary conditions
 * and a singular matrix. It uses a special approach to handle the singularity by fixing
 * the first element of the solution and recursively computing the rest.
 *
 * @param n The size of the matrix (number of rows/columns).
 * @param f The right-hand side vector.
 * @param x The solution vector. Output.
 *
 * @note
 * In the current version, p_const is forced to be 0.0, but actually, if you like,
 * you can set it to be x = x + p_const. This method handles the special case of
 * singular periodic systems by ensuring the solution satisfies the necessary constraints.
 */
void ChasingMethodBase::chasing_periodic_singular(int n, double* f, double* x)
{
    x[0] = 0.0; // In the current version, it is forced to be 0.0.
    x[1] = f[0] / n;
    for (int i = 1; i < n - 1; i++)
    {
        x[1] = x[1] - (n - 1 - i) * f[i] / n;
    }
    for (int i = 2; i < n; i++)
    {
        x[i] = f[i - 1] - x[i - 2] + 2.0 * x[i - 1];
    }
}
#include "gmres.h"

#include <cmath>
#include <iostream>
#include <vector>

// Solves the upper triangular system H(0:j,0:j) y = g(0:j)
// H is an (m+1)*m matrix, but only the top-left (j+1)x(j+1) block is used (upper triangular)
// j represents the dimension of the subproblem (starting from 0)
void back_substitution(const std::vector<double>& H, const std::vector<double>& g, std::vector<double>& y, int j)
{
    // j is the final subspace dimension; solve H(0:j,0:j) * y = g(0:j)
    // Note the storage format of H: H has size (m+1)×m
    // but only the top-left (j+1)×(j+1) block is used (i.e., H(i,j), i=0..j, j=0..j)
    // Row-major storage is used: the element H(i,j) is stored at H[i * m + j]
    for (int i = j; i >= 0; i--)
    {
        double sum = g[i];
        for (int k = i + 1; k <= j; k++)
        {
            sum -= H[i * (j + 1) + k] * y[k];
        }
        y[i] = sum / H[i * (j + 1) + i];
    }
}

// GMRES function
// Input:
// Afun: operator that performs multiplication with matrix A, assumed to be defined
// b: right-hand side vector
// m: dimension of the Krylov subspace
// tol: residual tolerance
// maxIter: maximum number of outer iterations
// x0: initial guess
// Output:
// x: solution vector
// resVec: history of residuals

field2 Afun(Shur_mat& S_1, Shur_mat& S_2, PoissonSolver2DInterface& solver_root, field2& x)
{
    field2 ft(x);
    ft = S_1 * x + S_2 * x;
    solver_root.solve(ft);
    return x - ft;
}

field2 Afun(Shur_mat& S_1, Shur_mat& S_2, Shur_mat& S_3, PoissonSolver2DInterface& solver_root, field2& x)
{
    field2 ft(x);
    ft = S_1 * x + S_2 * x + S_3 * x;
    solver_root.solve(ft);
    return x - ft;
}

field2
Afun(Shur_mat& S_1, Shur_mat& S_2, Shur_mat& S_3, Shur_mat& S_4, PoissonSolver2DInterface& solver_root, field2& x)
{
    field2 ft(x);
    ft = S_1 * x + S_2 * x + S_3 * x + S_4 * x;
    solver_root.solve(ft);
    return x - ft;
}

field2 Afun(std::vector<Shur_mat*>& S_params, PoissonSolver2DInterface& solver_root, field2& x)
{
    field2 ft(x);
    ft.clear(0.);
    for (auto& s : S_params)
    {
        ft = ft + *s * x;
    }
    // ft = S_1 * x + S_2 * x;
    solver_root.solve(ft);
    return x - ft;
}

field2 gmres(field2&                   b,
             field2&                   x0,
             Shur_mat&                 S_1,
             Shur_mat&                 S_2,
             PoissonSolver2DInterface& solver_root,
             int                       m,
             double                    tol,
             int                       maxIter,
             std::vector<double>&      resVec)
{
    // Initial residual r = b - A*x
    // std::vector<double> Ax(n,0.0);
    // Ax = Afun(x);

    // std::vector<double> r(n,0.0);

    field2 r(b);
    field2 x(x0);
    x           = x0;
    r           = b - Afun(S_1, S_2, solver_root, x);
    double beta = r.norm();
    if (beta < tol)
    {
        resVec.push_back(beta);
        return x;
    }

    field2 V[m + 1];
    for (int i = 0; i <= m; i++)
        V[i] = x;
    field2 w(x);

    for (int outer = 0; outer < maxIter; outer++)
    {
        // Allocate V and H
        // V: size n x (m+1)
        // std::vector<double> V(n*(m+1),0.0);
        // H: size (m+1) x m
        std::vector<double> H((m + 1) * m, 0.0);

        std::vector<double> cs(m, 0.0);
        std::vector<double> sn(m, 0.0);
        std::vector<double> g(m + 1, 0.0);

        // V(:,0) = r/beta
        V[0] = r * (1.0 / beta);
        g[0] = beta;

        for (int j = 0; j < m; j++)
        {
            w = Afun(S_1, S_2, solver_root, V[j]);

            {
                // Arnoldi Orthogonalization
                for (int i = 0; i <= j; i++)
                {
                    // H(i,j) = V(:,i)' * w
                    double hij   = V[i].dot(w);
                    H[i * m + j] = hij;
                    // w = w - H(i,j)*V(:,i)
                    w = w - V[i] * hij;
                }
                double h_j1j       = w.norm();
                H[(j + 1) * m + j] = h_j1j;
                if (h_j1j != 0)
                    V[j + 1] = w * (1.0 / h_j1j);
            }

            {
                // Apply the previous Givens rotation to the j-th column of H
                for (int i = 0; i < j; i++)
                {
                    double temp        = cs[i] * H[i * m + j] + sn[i] * H[(i + 1) * m + j];
                    H[(i + 1) * m + j] = -sn[i] * H[i * m + j] + cs[i] * H[(i + 1) * m + j];
                    H[i * m + j]       = temp;
                }
            }

            {
                // Construct new Givens rotation
                double h_jj          = H[j * m + j];
                double h_j1j_current = H[(j + 1) * m + j];
                double denom         = std::sqrt(h_jj * h_jj + h_j1j_current * h_j1j_current);
                cs[j]                = h_jj / denom;
                sn[j]                = h_j1j_current / denom;
            }

            {
                // Apply the new Givens rotation
                H[j * m + j]       = cs[j] * H[j * m + j] + sn[j] * H[(j + 1) * m + j];
                H[(j + 1) * m + j] = 0.0;
            }

            {
                // Update g vector
                double temp = cs[j] * g[j] + sn[j] * g[j + 1];
                g[j + 1]    = -sn[j] * g[j] + cs[j] * g[j + 1];
                g[j]        = temp;
            }

            double curRes = std::fabs(g[j + 1]);
            resVec.push_back(curRes);
            if (curRes < tol)
            {
                // Converged: solve upper triangular system for y
                std::vector<double> y(j + 1, 0.0);
                // Extract H(0:j,0:j)
                // Note that a j+1 x j+1 matrix is required
                // Reconstruct sub upper triangular matrix H(j+1)*(j+1)
                std::vector<double> H_sub((j + 1) * (j + 1), 0.0);
                for (int row = 0; row <= j; row++)
                {
                    for (int col = 0; col <= j; col++)
                    {
                        H_sub[row * (j + 1) + col] = H[row * m + col];
                    }
                }

                {
                    std::vector<double> g_sub(j + 1, 0.0);
                    for (int idx = 0; idx <= j; idx++)
                        g_sub[idx] = g[idx];

                    // Back substitution to solve y
                    // Note that H_sub is upper triangular
                    for (int row = j; row >= 0; row--)
                    {
                        double sum = g_sub[row];
                        for (int k = row + 1; k <= j; k++)
                        {
                            sum -= H_sub[row * (j + 1) + k] * y[k];
                        }
                        y[row] = sum / H_sub[row * (j + 1) + row];
                        // std::cout << "TP1-" << row << "-y=" << y[row] << std::endl;
                    }
                }

                {
                    // Update x = x + V(:,0:j)*y
                    for (int i = 0; i <= j; i++)
                        x = x + V[i] * y[i];
                }
                // std::cout << "=== TP1-" << "-x ==="<< std::endl;
                // x.print();
                // std::cout << "=============="<< std::endl;
                return x;
            }
        }

        // Inner loop ends, still not satisfied, solve H(0:m-1,0:m-1)*y=g(0:m-1)
        std::vector<double> y(m, 0.0);
        // Extract H(0:m,0:m)
        // Actually a m*m upper Hessenberg matrix, we need to solve H(0:m-1,0:m-1)*y = g(0:m-1)
        std::vector<double> H_sub(m * m, 0.0);
        for (int row = 0; row < m; row++)
        {
            for (int col = 0; col < m; col++)
            {
                H_sub[row * m + col] = H[row * m + col];
            }
        }
        std::vector<double> g_sub(m, 0.0);
        for (int idx = 0; idx < m; idx++)
            g_sub[idx] = g[idx];

        {
            // Solve y by back substitution
            for (int row = m - 1; row >= 0; row--)
            {
                double sum = g_sub[row];
                for (int k = row + 1; k < m; k++)
                {
                    sum -= H_sub[row * m + k] * y[k];
                }
                y[row] = sum / H_sub[row * m + row];
                // std::cout << "TP2-" << row << "-y=" << y[row] << std::endl;
            }
        }

        {
            // Update x = x + V(:,0:m-1)*y
            for (int i = 0; i < m; i++)
                x = x + V[i] * y[i];
        }
        // std::cout << "=== TP2-" << "-x ==="<< std::endl;
        // x.print();
        // std::cout << "=============="<< std::endl;

        {
            // Update residual
            r    = b - Afun(S_1, S_2, solver_root, x);
            beta = r.norm();
            resVec.push_back(beta);
        }

        if (beta < tol)
        {
            return x;
        }
    }
    return x;
}

field2 gmres(field2&                   b,
             field2&                   x0,
             Shur_mat&                 S_1,
             Shur_mat&                 S_2,
             Shur_mat&                 S_3,
             Shur_mat&                 S_4,
             PoissonSolver2DInterface& solver_root,
             int                       m,
             double                    tol,
             int                       maxIter,
             std::vector<double>&      resVec)
{
    field2 x(x0);
    double beta = 0.0;

    for (int outer = 0; outer < maxIter; outer++)
    {
        {
            field2 r = b - Afun(S_1, S_2, S_3, S_4, solver_root, x);
            beta     = r.norm();
            resVec.push_back(beta);
            if (beta < tol)
                return x;
        }

        std::vector<field2> V(m + 1, x0);        // Store basis vectors
        std::vector<double> H((m + 1) * m, 0.0); // Hessenberg matrix
        std::vector<double> cs(m, 0.0), sn(m, 0.0), g(m + 1, 0.0);

        // Initialize the first basis vector
        V[0] = (b - Afun(S_1, S_2, S_3, S_4, solver_root, x)) * (1.0 / beta);
        g[0] = beta;

        int j = 0;
        for (; j < m; j++)
        {
            {
                field2 w = Afun(S_1, S_2, S_3, S_4, solver_root, V[j]);

                // Orthogonalization
                for (int i = 0; i <= j; i++)
                {
                    H[i * m + j] = V[i].dot(w);
                    w            = w - V[i] * H[i * m + j];
                }

                const double h_j1j = w.norm();
                if (h_j1j < 1e-12)
                    break; // Early termination condition

                H[(j + 1) * m + j] = h_j1j;
                V[j + 1]           = w * (1.0 / h_j1j);
            }

            // Apply Givens rotation
            for (int i = 0; i < j; i++)
            {
                const double temp  = cs[i] * H[i * m + j] + sn[i] * H[(i + 1) * m + j];
                H[(i + 1) * m + j] = -sn[i] * H[i * m + j] + cs[i] * H[(i + 1) * m + j];
                H[i * m + j]       = temp;
            }

            // Compute new rotation
            const double h_jj  = H[j * m + j];
            const double h_j1j = H[(j + 1) * m + j];
            const double denom = std::hypot(h_jj, h_j1j);
            cs[j]              = h_jj / denom;
            sn[j]              = h_j1j / denom;

            // Update H and g
            H[j * m + j]       = cs[j] * H[j * m + j] + sn[j] * H[(j + 1) * m + j];
            H[(j + 1) * m + j] = 0.0;

            const double temp = cs[j] * g[j];
            g[j + 1]          = -sn[j] * g[j];
            g[j]              = temp;

            if (std::abs(g[j + 1]) < tol)
            {
                j++; // Include current j value
                break;
            }
        }

        // Back substitution to solve y
        std::vector<double> y(j);
        for (int i = j - 1; i >= 0; i--)
        {
            y[i] = g[i];
            for (int k = i + 1; k < j; k++)
                y[i] -= H[i * m + k] * y[k];
            y[i] /= H[i * m + i];
        }

        // Update solution
        for (int i = 0; i < j; i++)
            x = x + V[i] * y[i];
    }
    return x;
}

field2 gmres(field2&                   b,
             field2&                   x0,
             Shur_mat&                 S_1,
             Shur_mat&                 S_2,
             Shur_mat&                 S_3,
             PoissonSolver2DInterface& solver_root,
             int                       m,
             double                    tol,
             int                       maxIter,
             std::vector<double>&      resVec)
{
    // Initial residual r = b - A*x
    // std::vector<double> Ax(n,0.0);
    // Ax = Afun(x);

    // std::vector<double> r(n,0.0);

    field2 r(b);
    field2 x(x0);
    x           = x0;
    r           = b - Afun(S_1, S_2, S_3, solver_root, x);
    double beta = r.norm();
    if (beta < tol)
    {
        resVec.push_back(beta);
        return x;
    }

    field2 V[m + 1];
    for (int i = 0; i <= m; i++)
        V[i] = x;
    field2 w(x);

    for (int outer = 0; outer < maxIter; outer++)
    {
        // Allocate V and H
        // V: size n x (m+1)
        // std::vector<double> V(n*(m+1),0.0);
        // H: size (m+1) x m
        std::vector<double> H((m + 1) * m, 0.0);

        std::vector<double> cs(m, 0.0);
        std::vector<double> sn(m, 0.0);
        std::vector<double> g(m + 1, 0.0);

        // V(:,0) = r/beta
        V[0] = r * (1.0 / beta);
        g[0] = beta;

        for (int j = 0; j < m; j++)
        {
            w = Afun(S_1, S_2, S_3, solver_root, V[j]);

            {
                // Arnoldi Orthogonalization
                for (int i = 0; i <= j; i++)
                {
                    // H(i,j) = V(:,i)' * w
                    double hij   = V[i].dot(w);
                    H[i * m + j] = hij;
                    // w = w - H(i,j)*V(:,i)
                    w = w - V[i] * hij;
                }
                double h_j1j       = w.norm();
                H[(j + 1) * m + j] = h_j1j;
                if (h_j1j != 0)
                    V[j + 1] = w * (1.0 / h_j1j);
            }

            {
                // Apply the previous Givens rotation to the j-th column of H
                for (int i = 0; i < j; i++)
                {
                    double temp        = cs[i] * H[i * m + j] + sn[i] * H[(i + 1) * m + j];
                    H[(i + 1) * m + j] = -sn[i] * H[i * m + j] + cs[i] * H[(i + 1) * m + j];
                    H[i * m + j]       = temp;
                }
            }

            {
                // Construct new Givens rotation
                double h_jj          = H[j * m + j];
                double h_j1j_current = H[(j + 1) * m + j];
                double denom         = std::sqrt(h_jj * h_jj + h_j1j_current * h_j1j_current);
                cs[j]                = h_jj / denom;
                sn[j]                = h_j1j_current / denom;
            }

            {
                // Apply the new Givens rotation
                H[j * m + j]       = cs[j] * H[j * m + j] + sn[j] * H[(j + 1) * m + j];
                H[(j + 1) * m + j] = 0.0;
            }

            {
                // Update g vector
                double temp = cs[j] * g[j] + sn[j] * g[j + 1];
                g[j + 1]    = -sn[j] * g[j] + cs[j] * g[j + 1];
                g[j]        = temp;
            }

            double curRes = std::fabs(g[j + 1]);
            resVec.push_back(curRes);
            if (curRes < tol)
            {
                // Converged: solve upper triangular system for y
                std::vector<double> y(j + 1, 0.0);
                // Extract H(0:j,0:j)
                // Note that a j+1 x j+1 matrix is required
                // Reconstruct sub upper triangular matrix H(j+1)*(j+1)
                std::vector<double> H_sub((j + 1) * (j + 1), 0.0);
                for (int row = 0; row <= j; row++)
                {
                    for (int col = 0; col <= j; col++)
                    {
                        H_sub[row * (j + 1) + col] = H[row * m + col];
                    }
                }

                {
                    std::vector<double> g_sub(j + 1, 0.0);
                    for (int idx = 0; idx <= j; idx++)
                        g_sub[idx] = g[idx];
                    // Back substitution to solve y
                    // Note that H_sub is upper triangular
                    for (int row = j; row >= 0; row--)
                    {
                        double sum = g_sub[row];
                        for (int k = row + 1; k <= j; k++)
                        {
                            sum -= H_sub[row * (j + 1) + k] * y[k];
                        }
                        y[row] = sum / H_sub[row * (j + 1) + row];
                        // std::cout << "TP1-" << row << "-y=" << y[row] << std::endl;
                    }
                }

                {
                    // Update x = x + V(:,0:j)*y
                    for (int i = 0; i <= j; i++)
                        x = x + V[i] * y[i];
                }
                // std::cout << "=== TP1-" << "-x ==="<< std::endl;
                // x.print();
                // std::cout << "=============="<< std::endl;
                return x;
            }
        }

        // Inner loop ends, still not satisfied, solve H(0:m-1,0:m-1)*y=g(0:m-1)
        std::vector<double> y(m, 0.0);
        // Extract H(0:m,0:m)
        // Actually a m*m upper Hessenberg matrix, we need to solve H(0:m-1,0:m-1)*y = g(0:m-1)
        std::vector<double> H_sub(m * m, 0.0);
        for (int row = 0; row < m; row++)
        {
            for (int col = 0; col < m; col++)
            {
                H_sub[row * m + col] = H[row * m + col];
            }
        }
        std::vector<double> g_sub(m, 0.0);
        for (int idx = 0; idx < m; idx++)
            g_sub[idx] = g[idx];

        {
            // Solve y by back substitution
            for (int row = m - 1; row >= 0; row--)
            {
                double sum = g_sub[row];
                for (int k = row + 1; k < m; k++)
                {
                    sum -= H_sub[row * m + k] * y[k];
                }
                y[row] = sum / H_sub[row * m + row];
                // std::cout << "TP2-" << row << "-y=" << y[row] << std::endl;
            }
        }

        {
            // Update x = x + V(:,0:m-1)*y
            for (int i = 0; i < m; i++)
                x = x + V[i] * y[i];
        }
        // std::cout << "=== TP2-" << "-x ==="<< std::endl;
        // x.print();
        // std::cout << "=============="<< std::endl;

        {
            // Update residual
            r    = b - Afun(S_1, S_2, S_3, solver_root, x);
            beta = r.norm();
            resVec.push_back(beta);
        }

        if (beta < tol)
        {
            return x;
        }
    }
    return x;
}

field2 gmres(field2&                   b,
             field2&                   x0,
             std::vector<Shur_mat*>&   S_params,
             PoissonSolver2DInterface& solver_root,
             int                       m,
             double                    tol,
             int                       maxIter,
             std::vector<double>&      resVec)
{
    field2 x(x0);
    double beta = 0.0;

    for (int outer = 0; outer < maxIter; outer++)
    {
        {
            field2 r = b - Afun(S_params, solver_root, x);
            beta     = r.norm();
            resVec.push_back(beta);
            if (beta < tol)
                return x;
        }

        std::vector<field2> V(m + 1, x0);        // Store basis vectors
        std::vector<double> H((m + 1) * m, 0.0); // Hessenberg matrix
        std::vector<double> cs(m, 0.0), sn(m, 0.0), g(m + 1, 0.0);

        // Initialize the first basis vector
        V[0] = (b - Afun(S_params, solver_root, x)) * (1.0 / beta);
        g[0] = beta;

        int j = 0;
        for (; j < m; j++)
        {
            {
                field2 w = Afun(S_params, solver_root, V[j]);

                // Orthogonalization
                for (int i = 0; i <= j; i++)
                {
                    H[i * m + j] = V[i].dot(w);
                    w            = w - V[i] * H[i * m + j];
                }

                const double h_j1j = w.norm();
                if (h_j1j < 1e-12)
                    break; // Early termination condition

                H[(j + 1) * m + j] = h_j1j;
                V[j + 1]           = w * (1.0 / h_j1j);
            }

            // Apply Givens rotation
            for (int i = 0; i < j; i++)
            {
                const double temp  = cs[i] * H[i * m + j] + sn[i] * H[(i + 1) * m + j];
                H[(i + 1) * m + j] = -sn[i] * H[i * m + j] + cs[i] * H[(i + 1) * m + j];
                H[i * m + j]       = temp;
            }

            // Compute new rotation
            const double h_jj  = H[j * m + j];
            const double h_j1j = H[(j + 1) * m + j];
            const double denom = std::hypot(h_jj, h_j1j);
            cs[j]              = h_jj / denom;
            sn[j]              = h_j1j / denom;

            // Update H and g
            H[j * m + j]       = cs[j] * H[j * m + j] + sn[j] * H[(j + 1) * m + j];
            H[(j + 1) * m + j] = 0.0;

            const double temp = cs[j] * g[j];
            g[j + 1]          = -sn[j] * g[j];
            g[j]              = temp;

            if (std::abs(g[j + 1]) < tol)
            {
                j++; // Include current j value
                break;
            }
        }

        // Back substitution to solve y
        std::vector<double> y(j);
        for (int i = j - 1; i >= 0; i--)
        {
            y[i] = g[i];
            for (int k = i + 1; k < j; k++)
                y[i] -= H[i * m + k] * y[k];
            y[i] /= H[i * m + i];
        }

        // Update solution
        for (int i = 0; i < j; i++)
            x = x + V[i] * y[i];
    }
    return x;
}

// // -----------------
// Example

// // Example Afun: A is a diagonal matrix D(i,i)=i+1
// std::vector<double> Afun(const std::vector<double>& v) {
//     std::vector<double> w(v.size());
//     for (size_t i=0; i<v.size(); i++){
//         w[i] = (i+1)*v[i];
//     }
//     return w;
// }

// int main(){
//     size_t n = 5;
//     std::vector<double> b(n,1.0);
//     std::vector<double> x0(n,0.0);

//     int m = 5;
//     double tol = 1e-12;
//     int maxIter = 50;

//     std::vector<double> resVec;
//     std::vector<double> x = my_gmres(b, x0, m, tol, maxIter, resVec);

//     return 0;
// }

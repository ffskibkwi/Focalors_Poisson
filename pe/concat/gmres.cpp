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
field2 Afun(std::vector<SchurMat2D*>& S_params, DomainSolver2D& solver_root, field2& x)
{
    field2 ft(x);
    ft.clear(0.);
    for (auto& s : S_params)
        ft = ft + *s * x;
    solver_root.solve(ft);
    return x - ft;
}

field2 gmres(field2&                   b,
             field2&                   x0,
             std::vector<SchurMat2D*>&  S_params,
             DomainSolver2D& solver_root,
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

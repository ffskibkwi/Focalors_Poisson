#include "Shur_solver.hpp"

field2 Shur_solve_1(field2& f, Shur_mat& S, PoissonSolver2DInterface& solver_root, unsigned int max_it = 20)
{
    // Solve the equation: (A-S)x=f
    // 1:1 neibour
    field2 x_curr(f);
    field2 f_curr(f);

    solver_root.solve(f);
    x_curr = f;
    for (int it = 0; it < max_it; it++)
    {
        f_curr = S * x_curr;
        solver_root.solve(f_curr);
        x_curr = f + f_curr;
    }
    return x_curr;
}

field2
Shur_solve_2(field2& f, Shur_mat& S_1, Shur_mat& S_2, PoissonSolver2DInterface& solver_root, unsigned int max_it = 20)
{
    // Solve the equation: (A-S)x=f
    // 2:2 neibour
    field2 x_curr(f);
    field2 f_curr(f);

    solver_root.solve(f);
    x_curr = f;
    for (int it = 0; it < max_it; it++)
    {
        f_curr = S_1 * x_curr + S_2 * x_curr;
        solver_root.solve(f_curr);
        x_curr = f + f_curr;
    }
    return x_curr;
}

field2 Shur_solve_3(field2&                   f,
                    Shur_mat&                 S_1,
                    Shur_mat&                 S_2,
                    Shur_mat&                 S_3,
                    PoissonSolver2DInterface& solver_root,
                    unsigned int              max_it = 20)
{
    // Solve the equation: (A-S)x=f
    // 3:3 neibour
    field2 x_curr(f);
    field2 f_curr(f);

    solver_root.solve(f);
    x_curr = f;
    for (int it = 0; it < max_it; it++)
    {
        f_curr = S_1 * x_curr + S_2 * x_curr + S_3 * x_curr;
        solver_root.solve(f_curr);
        x_curr = f + f_curr;
    }
    return x_curr;
}

#pragma once

#include "Shur_mat.hpp"
#include "pe/poisson_old/poisson_solver_interface.h"

field2 Shur_solve_1(field2& f, Shur_mat& S, PoissonSolver2DInterface& solver_root, unsigned int max_it);
field2
Shur_solve_2(field2& f, Shur_mat& S_1, Shur_mat& S_2, PoissonSolver2DInterface& solver_root, unsigned int max_it);
field2 Shur_solve_3(field2&                   f,
                    Shur_mat&                 S_1,
                    Shur_mat&                 S_2,
                    Shur_mat&                 S_3,
                    PoissonSolver2DInterface& solver_root,
                    unsigned int              max_it);
#pragma once

#include "pch.h"

#include "Shur_mat.hpp"
#include "pe/poisson_old/poisson_solver_interface.h"

#include <vector>

field2 gmres(field2&                   b,
             field2&                   x,
             Shur_mat&                 S_1,
             Shur_mat&                 S_2,
             PoissonSolver2DInterface& solver_root,
             int                       m,
             double                    tol,
             int                       maxIter,
             std::vector<double>&      resVec);

field2 gmres(field2&                   b,
             field2&                   x0,
             Shur_mat&                 S_1,
             Shur_mat&                 S_2,
             Shur_mat&                 S_3,
             PoissonSolver2DInterface& solver_root,
             int                       m,
             double                    tol,
             int                       maxIter,
             std::vector<double>&      resVec);

field2 gmres(field2&                   b,
             field2&                   x,
             Shur_mat&                 S_1,
             Shur_mat&                 S_2,
             Shur_mat&                 S_3,
             Shur_mat&                 S_4,
             PoissonSolver2DInterface& solver_root,
             int                       m,
             double                    tol,
             int                       maxIter,
             std::vector<double>&      resVec);

field2 gmres(field2&                   b,
             field2&                   x,
             std::vector<Shur_mat*>&   S_params,
             PoissonSolver2DInterface& solver_root,
             int                       m,
             double                    tol,
             int                       maxIter,
             std::vector<double>&      resVec);
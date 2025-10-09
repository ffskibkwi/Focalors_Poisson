#pragma once

#include "pch.h"

#include "Schur_mat.h"
#include "pe/poisson/poisson_solver_interface.h"

#include <vector>

field2 gmres(field2&                   b,
             field2&                   x,
             std::vector<Schur_mat*>&   S_params,
             PoissonSolver2DInterface& solver_root,
             int                       m,
             double                    tol,
             int                       maxIter,
             std::vector<double>&      resVec);
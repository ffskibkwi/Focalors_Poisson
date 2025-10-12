#pragma once

#include "pch.h"

#include "Schur_mat.h"
// #include "pe/poisson_new/poisson_solver2d.h"
#include "domain_solver.h"

#include <vector>

field2 gmres(field2&                   b,
             field2&                   x0,
             std::vector<Schur_mat*>&  S_params,
             DomainSolver2D&           solver_root,
             int                       m,
             double                    tol,
             int                       maxIter,
             std::vector<double>&      resVec);
#pragma once

#include "base/pch.h"

#include "schur_mat2d.h"
// #include "pe/poisson/poisson_solver2d.h"
#include "domain_solver.h"

#include <vector>

field2 gmres(field2&                   b,
             field2&                   x0,
             std::vector<SchurMat2D*>&  S_params,
             DomainSolver2D&           solver_root,
             int                       m,
             double                    tol,
             int                       maxIter,
             std::vector<double>&      resVec);
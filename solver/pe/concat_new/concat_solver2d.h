#pragma once

#include "pch.h"
#include "core/boundary/boundary_type.h"
#include "core/domain/domain2d.h"
#include "core/domain/geometry2d.h"
#include "core/domain/variable.h"
#include "core/domain/geometry_tree.hpp"

#include "domain_solver.h"
#include "pe/poisson_new/poisson_solver2d.h"
#include "gmres_solver2d.h"

#include "Schur_mat.h"
#include <vector>
#include <unordered_map>

class ConcatSolver2D
{
    //Simple: Only for single main domain geometry
public:
    Variable& variable;
    
    ConcatSolver2D(Variable &in_variable);
    ~ConcatSolver2D();

    // void init();
    void solve();

private:
    std::unordered_map<Domain2DUniform*, field2*> temp_fields;
    
    std::unordered_map<Domain2DUniform*, DomainSolver2D*> solver_map;
    std::vector<Domain2DUniform*> solve_order;

    std::vector<double> resVec;

    void specify_solve_order();
    void construct_solver_map();
};

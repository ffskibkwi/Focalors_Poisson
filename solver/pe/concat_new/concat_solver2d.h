#pragma once

#include "pch.h"

#include "pe/poisson_new/poisson_solver2d.h"
#include "core/boundary/boundary_type.h"
#include "core/domain/domain2d.h"
#include "core/domain/geometry2d.h"
#include "core/domain/variable.h"
#include "core/domain/geometry_tree.hpp"
#include "Schur_mat.h"
#include <vector>
#include <stack>
#include <unordered_map>

class ConcatSolver2D
{
    //Simple: Only for single main domain geometry
public:
    Variable& variable;
    
    ConcatSolver2D(Variable &in_variable);
    ~ConcatSolver2D();

    void init();
    void prepare();
    void solve();

private:
    std::unordered_map<Domain2DUniform*, field2*> temp_fields;
    
    std::unordered_map<Domain2DUniform*, PoissonSolver2D*> pe_solver_map;
    std::unordered_map<Domain2DUniform*, DomainSolver2D*> solver_map;
    std::stack<Domain2DUniform*> solve_stack;
    
    // std::unordered_map<Domain2DUniform*, Schur_mat*> schur_mat_map;
    // std::vector<Schur_mat*> schur_mat_vec;
    // Domain2DUniform* geo_main_domain;

    std::vector<double> resVec;

    void construct_solve_queue();
};

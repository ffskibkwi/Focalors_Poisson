#pragma once

#include "base/pch.h"
#include "base/location_boundary.h"
#include "base/domain/domain2d.h"
#include "base/domain/geometry2d.h"
#include "base/domain/variable.h"
#include "base/domain/geometry_tree.hpp"

#include "domain_solver.h"
#include "pe/poisson/poisson_solver2d.h"
#include "gmres_solver2d.h"

#include "Schur_mat.h"
#include <vector>
#include <unordered_map>

class ConcatPoissonSolver2D
{
    //Simple: Only for single main domain geometry
public:
    Variable* variable = nullptr;
    
    ConcatPoissonSolver2D(Variable* in_variable, EnvironmentConfig* in_env_config = nullptr);
    ~ConcatPoissonSolver2D();

    // void init();
    void set_parameter(int in_m, double in_tol, int in_maxIter);
    void solve();

private:
    std::unordered_map<Domain2DUniform*, field2*> temp_fields;
    
    std::unordered_map<Domain2DUniform*, DomainSolver2D*> solver_map;
    std::vector<Domain2DUniform*> solve_order;

    std::vector<double> resVec;

    int    m       = 20;
    double tol     = 1e-8;
    int    maxIter = 100;

    void specify_solve_order();
    void construct_solver_map();

    Domain2DUniform* tree_root = nullptr;
    std::unordered_map<Domain2DUniform*, field2*> field_map;
    std::unordered_map<Domain2DUniform*, std::unordered_map<LocationType, Domain2DUniform*>> tree_map;
    std::unordered_map<Domain2DUniform*, std::pair<LocationType, Domain2DUniform*>> parent_map;

    EnvironmentConfig* env_config;
    bool showGmresRes = false;
};

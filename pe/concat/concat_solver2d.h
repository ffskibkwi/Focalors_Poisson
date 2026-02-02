#pragma once

#include "base/pch.h"

#include "base/domain/domain2d.h"
#include "base/domain/geometry2d.h"
#include "base/domain/geometry_tree.hpp"
#include "base/domain/variable2d.h"
#include "base/location_boundary.h"
#include "domain_solver.h"
#include "io/config.h"
#include "schur_mat2d.h"

#include <chrono>
#include <unordered_map>
#include <vector>

class ConcatPoissonSolver2D
{
public:
    Variable2D* variable = nullptr;

    ConcatPoissonSolver2D() {}
    ConcatPoissonSolver2D(Variable2D* in_variable, EnvironmentConfig* in_env_config = nullptr);
    ~ConcatPoissonSolver2D();

    void set_parameter(int in_m, double in_tol, int in_maxIter);
    void solve();

protected:
    void init_before_constructing_solver(Variable2D* _variable);

    void construct_solver_map_at_domain(Domain2DUniform* domain);
    void construct_solver_map();

    void boundary_assembly();

    std::unordered_map<Domain2DUniform*, field2*> temp_fields;

    std::unordered_map<Domain2DUniform*, DomainSolver2D*> solver_map;
    std::vector<std::vector<Domain2DUniform*>>            hierarchical_solve_levels;

    std::vector<double> resVec;

    int    m       = 20;
    double tol     = 1e-8;
    int    maxIter = 100;

    Domain2DUniform*                                                                         tree_root = nullptr;
    std::unordered_map<Domain2DUniform*, field2*>                                            field_map;
    std::unordered_map<Domain2DUniform*, std::unordered_map<LocationType, Domain2DUniform*>> tree_map;
    std::unordered_map<Domain2DUniform*, std::pair<LocationType, Domain2DUniform*>>          parent_map;

    bool                                      showGmresRes = false;
    EnvironmentConfig*                        env_config;
    bool                                      track_time = false;
    std::chrono::duration<double, std::milli> schur_total; // should be clear at the beginning of construct_solver_map
};

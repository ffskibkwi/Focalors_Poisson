#pragma once

#include "base/pch.h"

#include "base/domain/domain3d.h"
#include "base/domain/geometry3d.h"
#include "base/domain/geometry_tree.hpp"
#include "base/domain/variable2d.h"
#include "base/domain/variable3d.h"
#include "base/location_boundary.h"
#include "domain_solver.h"
#include "io/config.h"

#include "schur_mat3d.h"
#include <unordered_map>
#include <vector>

class ConcatPoissonSolver3D
{
public:
    Variable3D* variable = nullptr;

    ConcatPoissonSolver3D(Variable3D* in_variable, EnvironmentConfig* in_env_config = nullptr);
    ~ConcatPoissonSolver3D();

    void set_parameter(int in_m, double in_tol, int in_maxIter);
    void solve();

private:
    void construct_solver_map_at_domain(Domain3DUniform* domain);
    void construct_solver_map();

    void boundary_assembly();

    std::unordered_map<Domain3DUniform*, field3*> temp_fields;

    std::unordered_map<Domain3DUniform*, DomainSolver3D*> solver_map;
    std::vector<std::vector<Domain3DUniform*>>            hierarchical_solve_levels;

    std::vector<double> resVec;

    int    m       = 20;
    double tol     = 1e-8;
    int    maxIter = 100;

    Domain3DUniform*                                                                         tree_root = nullptr;
    std::unordered_map<Domain3DUniform*, field3*>                                            field_map;
    std::unordered_map<Domain3DUniform*, std::unordered_map<LocationType, Domain3DUniform*>> tree_map;
    std::unordered_map<Domain3DUniform*, std::pair<LocationType, Domain3DUniform*>>          parent_map;

    EnvironmentConfig* env_config;
    bool               showGmresRes = false;
    bool               track_time   = false;
};
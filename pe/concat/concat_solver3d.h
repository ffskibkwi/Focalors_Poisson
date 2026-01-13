#pragma once

#include "base/pch.h"
#include "base/location_boundary.h"
#include "base/domain/domain3d.h"
#include "base/domain/geometry3d.h"
#include "base/domain/variable.h"
#include "base/domain/geometry_tree.hpp"

#include "domain_solver.h"
#include "pe/poisson/poisson_solver3d.h"
#include "gmres_solver3d.h"

#include "Schur_mat3d.h"
#include <vector>
#include <unordered_map>

class ConcatPoissonSolver3D
{
    //Simple: Only for single main domain geometry
public:
    Variable* variable = nullptr;
    
    ConcatPoissonSolver3D(Variable* in_variable, EnvironmentConfig* in_env_config = nullptr);
    ~ConcatPoissonSolver3D();

    void set_parameter(int in_m, double in_tol, int in_maxIter);
    void solve();

private:
    std::unordered_map<Domain3DUniform*, field3*> temp_fields;
    
    std::unordered_map<Domain3DUniform*, DomainSolver3D*> solver_map;
    std::vector<Domain3DUniform*> solve_order;

    std::vector<double> resVec;

    int    m       = 20;
    double tol     = 1e-8;
    int    maxIter = 100;

    void specify_solve_order();
    void construct_solver_map();

    Domain3DUniform* tree_root = nullptr;
    std::unordered_map<Domain3DUniform*, field3*> field_map;
    std::unordered_map<Domain3DUniform*, std::unordered_map<LocationType, Domain3DUniform*>> tree_map;
    std::unordered_map<Domain3DUniform*, std::pair<LocationType, Domain3DUniform*>> parent_map;

    EnvironmentConfig* env_config;
    bool showGmresRes = false;

    void set_boundary();
    
    // 辅助函数：实现 field3 的 bond_add 功能
    void bond_add_3d(field3& target, LocationType location, double k, const field3& source);
};
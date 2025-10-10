#pragma once

#include "pch.h"

#include "pe/poisson_new/poisson_solver2d.h"
#include "core/boundary/boundary_type.h"
#include "core/domain/domain2d.h"
#include "core/domain/geometry2d.h"
#include "core/domain/variable.h"
#include "Schur_mat.h"
#include <vector>
#include <unordered_map>

class ConcatSolver2D_Simple
{
    //Simple: Only for single main domain geometry
public:
    Variable& variable;
    
    ConcatSolver2D_Simple(Variable &in_variable);
    ~ConcatSolver2D_Simple();

    void init();
    void solve();

private:
    std::unordered_map<Domain2DUniform*, PoissonSolver2D*> pe_solver_map;
    std::unordered_map<Domain2DUniform*, field2*> temp_fields;

    std::unordered_map<Domain2DUniform*, Schur_mat*> schur_mat_map;
    std::vector<Schur_mat*> schur_mat_vec;
    
    Domain2DUniform* geo_main_domain;

    std::unordered_map<LocationType, Domain2DUniform*> geo_main_neighbour_domains;  //Only for single main domain geometry

    std::vector<double> resVec;
};





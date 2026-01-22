#pragma once

#include "base/pch.h"
#include "base/location_boundary.h"
#include "pe/poisson_old/poisson_solver_interface.h"

class SchurMat2D
{
protected:
    int      cosize_n;
    int      root_nx, root_ny, branch_nx, branch_ny;
    double** value;

public:
    LocationType direction;

    SchurMat2D(const field2& root, const field2& branch, LocationType dir);
    ~SchurMat2D();

    void construct(PoissonSolver2DInterface& branch_solver);
    field2 operator*(const field2& root);

    void print();
};



#pragma once
#include "core/base/field2.h"
#include "core/base/field3.h"

// This provides base class for PoissonSolver/GMRESSolver
// This is used for solver queue in concat solver
// P.S. Although it named DomainSolver, it do not need to link with a domain

class DomainSolver2D
{
public:
    DomainSolver2D() {};

    ~DomainSolver2D();

    // virtual void init() = 0;
    virtual void solve(field2& f) = 0;
};

class DomainSolver3D
{
public:
    DomainSolver3D() {};

    ~DomainSolver3D();

    // virtual void init() = 0;
    virtual void solve(field3& f) = 0;
};
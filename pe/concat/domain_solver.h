#pragma once
#include "base/field/field2.h"
#include "base/field/field3.h"

// This provides base class for PoissonSolver/GMRESSolver
// This is used for solver queue in concat solver
// P.S. Although it named DomainSolver, it do not need to link with a domain

class DomainSolver2D
{
public:
    DomainSolver2D() {};
    virtual ~DomainSolver2D() {};
    virtual void solve(field2& f) = 0;
};

class DomainSolver3D
{
public:
    DomainSolver3D() {};
    virtual ~DomainSolver3D() {};
    virtual void solve(field3& f) = 0;
};
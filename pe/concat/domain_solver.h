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

    virtual ~DomainSolver2D() = default;

    // virtual void init() = 0;
    virtual void solve(field2& f) = 0;

    // 集体求解（约定：通信器 root 持有完整 b，其他进程仅参与通信）
    virtual void solve_collective_root_owned(field2& f) { solve(f); }

    // 是否为内部通信器 root（串行/无通信器默认 true）
    virtual bool is_comm_root() const { return true; }
};

class DomainSolver3D
{
public:
    DomainSolver3D() {};

    virtual ~DomainSolver3D() = default;

    // virtual void init() = 0;
    virtual void solve(field3& f) = 0;
};
#pragma once

#include "base/pch.h"

class DomainSolver3D;

class SchurMat3D
{
protected:
    int    bnx; // branch domain nx
    int    bny; // branch domain ny
    int    bnz; // branch domain nz
    int    cn;  // concat interface total n
    field2 value;

public:
    SchurMat3D(const Domain3DUniform* domain, const int _cn)
        : bnx(domain->get_nx())
        , bny(domain->get_ny())
        , bnz(domain->get_nz())
        , cn(_cn)
        , value(_cn, _cn)
    {}
    int  get_size() const { return cn; }
    void print() { value.print(); }

    virtual void   construct(DomainSolver3D* branch_solver) = 0;
    virtual field3 operator*(const field3& root)            = 0;

    void set_name(const std::string& name) { value.set_name(name); }
};

// Left face: x=0 plane, size = ny * nz
class SchurMat3D_left : public SchurMat3D
{
public:
    SchurMat3D_left(const Domain3DUniform* domain)
        : SchurMat3D(domain, domain->get_ny() * domain->get_nz())
    {}
    void   construct(DomainSolver3D* branch_solver) override;
    field3 operator*(const field3& root) override;
};

// Right face: x=nx-1 plane, size = ny * nz
class SchurMat3D_right : public SchurMat3D
{
public:
    SchurMat3D_right(const Domain3DUniform* domain)
        : SchurMat3D(domain, domain->get_ny() * domain->get_nz())
    {}
    void   construct(DomainSolver3D* branch_solver) override;
    field3 operator*(const field3& root) override;
};

// Front face: y=0 plane, size = nx * nz
class SchurMat3D_front : public SchurMat3D
{
public:
    SchurMat3D_front(const Domain3DUniform* domain)
        : SchurMat3D(domain, domain->get_nx() * domain->get_nz())
    {}
    void   construct(DomainSolver3D* branch_solver) override;
    field3 operator*(const field3& root) override;
};

// Back face: y=ny-1 plane, size = nx * nz
class SchurMat3D_back : public SchurMat3D
{
public:
    SchurMat3D_back(const Domain3DUniform* domain)
        : SchurMat3D(domain, domain->get_nx() * domain->get_nz())
    {}
    void   construct(DomainSolver3D* branch_solver) override;
    field3 operator*(const field3& root) override;
};

// Down face: z=0 plane, size = nx * ny
class SchurMat3D_down : public SchurMat3D
{
public:
    SchurMat3D_down(const Domain3DUniform* domain)
        : SchurMat3D(domain, domain->get_nx() * domain->get_ny())
    {}
    void   construct(DomainSolver3D* branch_solver) override;
    field3 operator*(const field3& root) override;
};

// Up face: z=nz-1 plane, size = nx * ny
class SchurMat3D_up : public SchurMat3D
{
public:
    SchurMat3D_up(const Domain3DUniform* domain)
        : SchurMat3D(domain, domain->get_nx() * domain->get_ny())
    {}
    void   construct(DomainSolver3D* branch_solver) override;
    field3 operator*(const field3& root) override;
};
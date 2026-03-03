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

// XNegative face: x=0 plane, size = ny * nz
class SchurMat3D_xneg : public SchurMat3D
{
public:
    SchurMat3D_xneg(const Domain3DUniform* domain)
        : SchurMat3D(domain, domain->get_ny() * domain->get_nz())
    {}
    void   construct(DomainSolver3D* branch_solver) override;
    field3 operator*(const field3& root) override;
};

// XPositive face: x=nx-1 plane, size = ny * nz
class SchurMat3D_xpos : public SchurMat3D
{
public:
    SchurMat3D_xpos(const Domain3DUniform* domain)
        : SchurMat3D(domain, domain->get_ny() * domain->get_nz())
    {}
    void   construct(DomainSolver3D* branch_solver) override;
    field3 operator*(const field3& root) override;
};

// YNegative face: y=0 plane, size = nx * nz
class SchurMat3D_yneg : public SchurMat3D
{
public:
    SchurMat3D_yneg(const Domain3DUniform* domain)
        : SchurMat3D(domain, domain->get_nx() * domain->get_nz())
    {}
    void   construct(DomainSolver3D* branch_solver) override;
    field3 operator*(const field3& root) override;
};

// YPositive face: y=ny-1 plane, size = nx * nz
class SchurMat3D_ypos : public SchurMat3D
{
public:
    SchurMat3D_ypos(const Domain3DUniform* domain)
        : SchurMat3D(domain, domain->get_nx() * domain->get_nz())
    {}
    void   construct(DomainSolver3D* branch_solver) override;
    field3 operator*(const field3& root) override;
};

// ZNegative face: z=0 plane, size = nx * ny
class SchurMat3D_zneg : public SchurMat3D
{
public:
    SchurMat3D_zneg(const Domain3DUniform* domain)
        : SchurMat3D(domain, domain->get_nx() * domain->get_ny())
    {}
    void   construct(DomainSolver3D* branch_solver) override;
    field3 operator*(const field3& root) override;
};

// ZPositive face: z=nz-1 plane, size = nx * ny
class SchurMat3D_zpos : public SchurMat3D
{
public:
    SchurMat3D_zpos(const Domain3DUniform* domain)
        : SchurMat3D(domain, domain->get_nx() * domain->get_ny())
    {}
    void   construct(DomainSolver3D* branch_solver) override;
    field3 operator*(const field3& root) override;
};
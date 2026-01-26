#pragma once

#include "base/pch.h"

class DomainSolver2D;

class SchurMat2D
{
protected:
    int    bnx; // branch domain nx
    int    bny; // branch domain ny
    int    cn;  // concat interface total n
    field2 value;

public:
    SchurMat2D(const Domain2DUniform* domain, const int _cn)
        : bnx(domain->get_nx())
        , bny(domain->get_ny())
        , cn(_cn)
        , value(_cn, _cn)
    {}
    int  get_size() const { return cn; }
    void print() { value.print(); }

    virtual void   construct(DomainSolver2D* branch_solver) = 0;
    virtual field2 operator*(const field2& root)            = 0;

    void set_name(const std::string& name) { value.set_name(name); }
    void write_csv(const std::string& directory);
};

class SchurMat2D_left : public SchurMat2D
{
public:
    SchurMat2D_left(const Domain2DUniform* domain)
        : SchurMat2D(domain, domain->get_ny())
    {}
    void   construct(DomainSolver2D* branch_solver) override;
    field2 operator*(const field2& root) override;
};

class SchurMat2D_right : public SchurMat2D
{
public:
    SchurMat2D_right(const Domain2DUniform* domain)
        : SchurMat2D(domain, domain->get_ny())
    {}
    void   construct(DomainSolver2D* branch_solver) override;
    field2 operator*(const field2& root) override;
};

class SchurMat2D_up : public SchurMat2D
{
public:
    SchurMat2D_up(const Domain2DUniform* domain)
        : SchurMat2D(domain, domain->get_nx())
    {}
    void   construct(DomainSolver2D* branch_solver) override;
    field2 operator*(const field2& root) override;
};

class SchurMat2D_down : public SchurMat2D
{
public:
    SchurMat2D_down(const Domain2DUniform* domain)
        : SchurMat2D(domain, domain->get_nx())
    {}
    void   construct(DomainSolver2D* branch_solver) override;
    field2 operator*(const field2& root) override;
};

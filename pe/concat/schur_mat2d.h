#pragma once

#include "base/pch.h"

#include "base/domain/domain2d.h"
#include "domain_solver.h"


class SchurMat2D
{
protected:
    int      cosize_n;
    int      root_nx, root_ny, branch_nx, branch_ny;
    double** value;

public:
    SchurMat2D(const field2& root, const field2& branch, const int _cosize_n)
        : cosize_n(_cosize_n)
        , root_nx(root.get_nx())
        , root_ny(root.get_ny())
        , branch_nx(branch.get_nx())
        , branch_ny(branch.get_ny())
    {
        value = new double*[cosize_n];
        for (int i = 0; i < cosize_n; i++)
            value[i] = new double[cosize_n];
    }
    SchurMat2D(const Domain2DUniform& root, const Domain2DUniform& branch, const int _cosize_n)
        : cosize_n(_cosize_n)
        , root_nx(root.nx)
        , root_ny(root.ny)
        , branch_nx(branch.nx)
        , branch_ny(branch.ny)
    {
        value = new double*[cosize_n];
        for (int i = 0; i < cosize_n; i++)
            value[i] = new double[cosize_n];
    }
    ~SchurMat2D()
    {
        for (int i = 0; i < cosize_n; i++)
            delete[] value[i];
        delete[] value;
    }
    int  get_size() const { return cosize_n; }
    void store_rowmajor(double* buf) const
    {
        for (int i = 0; i < cosize_n; ++i)
            for (int j = 0; j < cosize_n; ++j)
                buf[i * cosize_n + j] = value[i][j];
    }
    void load_rowmajor(const double* buf)
    {
        for (int i = 0; i < cosize_n; ++i)
            for (int j = 0; j < cosize_n; ++j)
                value[i][j] = buf[i * cosize_n + j];
    }
    void print()
    {
        for (int i = 0; i < cosize_n; i++)
        {
            for (int j = 0; j < cosize_n; j++)
                std::cout << value[i][j] << std::endl;
        }
    }

    virtual void   construct(DomainSolver2D* branch_solver) = 0;
    virtual field2 operator*(const field2& root)            = 0;

    void set_name(const std::string& in_name) { name = in_name; }
    void dump_to_csv(const std::string& directory);

protected:
    std::string name = "SchurMat";
};

class SchurMat2D_left : public SchurMat2D
{
public:
    SchurMat2D_left(const Domain2DUniform& root, const Domain2DUniform& branch)
        : SchurMat2D(root, branch, root.ny)
    {}
    void   construct(DomainSolver2D* branch_solver) override;
    field2 operator*(const field2& root) override;
};

class SchurMat2D_right : public SchurMat2D
{
public:
    SchurMat2D_right(const Domain2DUniform& root, const Domain2DUniform& branch)
        : SchurMat2D(root, branch, root.ny)
    {}
    void   construct(DomainSolver2D* branch_solver) override;
    field2 operator*(const field2& root) override;
};

class SchurMat2D_up : public SchurMat2D
{
public:
    SchurMat2D_up(const Domain2DUniform& root, const Domain2DUniform& branch)
        : SchurMat2D(root, branch, root.nx)
    {}
    void   construct(DomainSolver2D* branch_solver) override;
    field2 operator*(const field2& root) override;
};

class SchurMat2D_down : public SchurMat2D
{
public:
    SchurMat2D_down(const Domain2DUniform& root, const Domain2DUniform& branch)
        : SchurMat2D(root, branch, root.nx)
    {}
    void   construct(DomainSolver2D* branch_solver) override;
    field2 operator*(const field2& root) override;
};

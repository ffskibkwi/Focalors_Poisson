#pragma once

#include "pch.h"

#include "pe/poisson/poisson_solver_interface.h"

class Shur_mat
{
protected:
    int      cosize_n;
    int      root_nx, root_ny, branch_nx, branch_ny;
    double** value;

public:
    Shur_mat(const field2& root, const field2& branch, const int _cosize_n)
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
    ~Shur_mat()
    {
        for (int i = 0; i < cosize_n; i++)
            delete[] value[i];
        delete[] value;
    }
    void print()
    {
        for (int i = 0; i < cosize_n; i++)
        {
            for (int j = 0; j < cosize_n; j++)
                std::cout << value[i][j] << std::endl;
        }
    }

    virtual void   construct(PoissonSolver2DInterface& branch_solver) = 0;
    virtual field2 operator*(const field2& root)                      = 0;
};

class Shur_mat_left : public Shur_mat
{
public:
    Shur_mat_left(const field2& root, const field2& branch)
        : Shur_mat(root, branch, root.get_ny())
    {}
    void   construct(PoissonSolver2DInterface& branch_solver) override;
    field2 operator*(const field2& root) override;
};

class Shur_mat_right : public Shur_mat
{
public:
    Shur_mat_right(const field2& root, const field2& branch)
        : Shur_mat(root, branch, root.get_ny())
    {}
    void   construct(PoissonSolver2DInterface& branch_solver) override;
    field2 operator*(const field2& root) override;
};

class Shur_mat_up : public Shur_mat
{
public:
    Shur_mat_up(const field2& root, const field2& branch)
        : Shur_mat(root, branch, root.get_nx())
    {}
    void   construct(PoissonSolver2DInterface& branch_solver) override;
    field2 operator*(const field2& root) override;
};

class Shur_mat_down : public Shur_mat
{
public:
    Shur_mat_down(const field2& root, const field2& branch)
        : Shur_mat(root, branch, root.get_nx())
    {}
    void   construct(PoissonSolver2DInterface& branch_solver) override;
    field2 operator*(const field2& root) override;
};

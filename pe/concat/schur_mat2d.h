#pragma once

#include "base/pch.h"

class DomainSolver2D;

class SchurMat2D
{
protected:
    int      branch_nx, branch_ny;
    int      concat_n;
    double** value;

public:
    SchurMat2D(const int _branch_nx, const int _branch_ny, const int _concat_n)
        : branch_nx(_branch_nx)
        , branch_ny(_branch_ny)
        , concat_n(_concat_n)
    {
        value = new double*[concat_n];
        for (int i = 0; i < concat_n; i++)
            value[i] = new double[concat_n];
    }
    ~SchurMat2D()
    {
        for (int i = 0; i < concat_n; i++)
            delete[] value[i];
        delete[] value;
    }
    int  get_size() const { return concat_n; }
    void print()
    {
        for (int i = 0; i < concat_n; i++)
            for (int j = 0; j < concat_n; j++)
                std::cout << value[i][j] << std::endl;
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
    SchurMat2D_left(const int _branch_nx, const int _branch_ny)
        : SchurMat2D(_branch_nx, _branch_ny, _branch_ny)
    {}
    void   construct(DomainSolver2D* branch_solver) override;
    field2 operator*(const field2& root) override;
};

class SchurMat2D_right : public SchurMat2D
{
public:
    SchurMat2D_right(const int _branch_nx, const int _branch_ny)
        : SchurMat2D(_branch_nx, _branch_ny, _branch_ny)
    {}
    void   construct(DomainSolver2D* branch_solver) override;
    field2 operator*(const field2& root) override;
};

class SchurMat2D_up : public SchurMat2D
{
public:
    SchurMat2D_up(const int _branch_nx, const int _branch_ny)
        : SchurMat2D(_branch_nx, _branch_ny, _branch_nx)
    {}
    void   construct(DomainSolver2D* branch_solver) override;
    field2 operator*(const field2& root) override;
};

class SchurMat2D_down : public SchurMat2D
{
public:
    SchurMat2D_down(const int _branch_nx, const int _branch_ny)
        : SchurMat2D(_branch_nx, _branch_ny, _branch_nx)
    {}
    void   construct(DomainSolver2D* branch_solver) override;
    field2 operator*(const field2& root) override;
};

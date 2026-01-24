#pragma once

#include "base/pch.h"

class DomainSolver2D;

class SchurMat2D
{
protected:
    int      bnx; // branch domain nx
    int      bny; // branch domain ny
    int      cn;  // concat interface total n
    double** value;

public:
    SchurMat2D(const int _bnx, const int _bny, const int _cn)
        : bnx(_bnx)
        , bny(_bny)
        , cn(_cn)
    {
        value = new double*[cn];
        for (int i = 0; i < cn; i++)
            value[i] = new double[cn];
    }
    ~SchurMat2D()
    {
        for (int i = 0; i < cn; i++)
            delete[] value[i];
        delete[] value;
    }
    int  get_size() const { return cn; }
    void print()
    {
        for (int i = 0; i < cn; i++)
            for (int j = 0; j < cn; j++)
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
    SchurMat2D_left(const int _bnx, const int _bny)
        : SchurMat2D(_bnx, _bny, _bny)
    {}
    void   construct(DomainSolver2D* branch_solver) override;
    field2 operator*(const field2& root) override;
};

class SchurMat2D_right : public SchurMat2D
{
public:
    SchurMat2D_right(const int _bnx, const int _bny)
        : SchurMat2D(_bnx, _bny, _bny)
    {}
    void   construct(DomainSolver2D* branch_solver) override;
    field2 operator*(const field2& root) override;
};

class SchurMat2D_up : public SchurMat2D
{
public:
    SchurMat2D_up(const int _bnx, const int _bny)
        : SchurMat2D(_bnx, _bny, _bnx)
    {}
    void   construct(DomainSolver2D* branch_solver) override;
    field2 operator*(const field2& root) override;
};

class SchurMat2D_down : public SchurMat2D
{
public:
    SchurMat2D_down(const int _bnx, const int _bny)
        : SchurMat2D(_bnx, _bny, _bnx)
    {}
    void   construct(DomainSolver2D* branch_solver) override;
    field2 operator*(const field2& root) override;
};

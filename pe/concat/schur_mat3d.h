#pragma once

#include "base/pch.h"

class DomainSolver3D;

class SchurMat3D
{
protected:
    int      bnx; // branch domain nx
    int      bny; // branch domain ny
    int      bnz; // branch domain nz
    int      cn;  // concat interface total n
    double** value;

public:
    SchurMat3D(const int _bnx, const int _bny, const int _bnz, const int _cn)
        : bnx(_bnx)
        , bny(_bny)
        , bnz(_bnz)
        , cn(_cn)
    {
        value = new double*[cn];
        for (int i = 0; i < cn; i++)
            value[i] = new double[cn];
    }
    ~SchurMat3D()
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

    virtual void   construct(DomainSolver3D* branch_solver) = 0;
    virtual field3 operator*(const field3& root)            = 0;
};

// Left face: x=0 plane, size = ny * nz
class SchurMat3D_left : public SchurMat3D
{
public:
    SchurMat3D_left(const int _bnx, const int _bny, const int _bnz)
        : SchurMat3D(_bnx, _bny, _bnz, _bny * _bnz)
    {}
    void   construct(DomainSolver3D* branch_solver) override;
    field3 operator*(const field3& root) override;
};

// Right face: x=nx-1 plane, size = ny * nz
class SchurMat3D_right : public SchurMat3D
{
public:
    SchurMat3D_right(const int _bnx, const int _bny, const int _bnz)
        : SchurMat3D(_bnx, _bny, _bnz, _bny * _bnz)
    {}
    void   construct(DomainSolver3D* branch_solver) override;
    field3 operator*(const field3& root) override;
};

// Front face: y=0 plane, size = nx * nz
class SchurMat3D_front : public SchurMat3D
{
public:
    SchurMat3D_front(const int _bnx, const int _bny, const int _bnz)
        : SchurMat3D(_bnx, _bny, _bnz, _bnx * _bnz)
    {}
    void   construct(DomainSolver3D* branch_solver) override;
    field3 operator*(const field3& root) override;
};

// Back face: y=ny-1 plane, size = nx * nz
class SchurMat3D_back : public SchurMat3D
{
public:
    SchurMat3D_back(const int _bnx, const int _bny, const int _bnz)
        : SchurMat3D(_bnx, _bny, _bnz, _bnx * _bnz)
    {}
    void   construct(DomainSolver3D* branch_solver) override;
    field3 operator*(const field3& root) override;
};

// Down face: z=0 plane, size = nx * ny
class SchurMat3D_down : public SchurMat3D
{
public:
    SchurMat3D_down(const int _bnx, const int _bny, const int _bnz)
        : SchurMat3D(_bnx, _bny, _bnz, _bnx * _bny)
    {}
    void   construct(DomainSolver3D* branch_solver) override;
    field3 operator*(const field3& root) override;
};

// Up face: z=nz-1 plane, size = nx * ny
class SchurMat3D_up : public SchurMat3D
{
public:
    SchurMat3D_up(const int _bnx, const int _bny, const int _bnz)
        : SchurMat3D(_bnx, _bny, _bnz, _bnx * _bny)
    {}
    void   construct(DomainSolver3D* branch_solver) override;
    field3 operator*(const field3& root) override;
};
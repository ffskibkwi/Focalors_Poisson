#pragma once

#include "base/pch.h"

class DomainSolver3D;

class SchurMat3D
{
protected:
    int      root_nx, root_ny, root_nz, branch_nx, branch_ny, branch_nz;
    int      concat_nx;
    double** value;

public:
    SchurMat3D(const int _root_nx,
               const int _root_ny,
               const int _root_nz,
               const int _branch_nx,
               const int _branch_ny,
               const int _branch_nz,
               const int _concat_nx)
        : root_nx(_root_nx)
        , root_ny(_root_ny)
        , root_nz(_root_nz)
        , branch_nx(_branch_nx)
        , branch_ny(_branch_ny)
        , branch_nz(_branch_nz)
        , concat_nx(_concat_nx)
    {
        value = new double*[concat_nx];
        for (int i = 0; i < concat_nx; i++)
            value[i] = new double[concat_nx];
    }
    ~SchurMat3D()
    {
        for (int i = 0; i < concat_nx; i++)
            delete[] value[i];
        delete[] value;
    }
    int  get_size() const { return concat_nx; }
    void print()
    {
        for (int i = 0; i < concat_nx; i++)
        {
            for (int j = 0; j < concat_nx; j++)
                std::cout << value[i][j] << std::endl;
        }
    }

    virtual void   construct(DomainSolver3D* branch_solver) = 0;
    virtual field3 operator*(const field3& root)            = 0;
};

// Left face: x=0 plane, size = ny * nz
class SchurMat3D_left : public SchurMat3D
{
public:
    SchurMat3D_left(const int _root_nx,
                    const int _root_ny,
                    const int _root_nz,
                    const int _branch_nx,
                    const int _branch_ny,
                    const int _branch_nz)
        : SchurMat3D(_root_nx, _root_ny, _root_nz, _branch_nx, _branch_ny, _branch_nz, _root_ny * _root_nz)
    {}
    void   construct(DomainSolver3D* branch_solver) override;
    field3 operator*(const field3& root) override;
};

// Right face: x=nx-1 plane, size = ny * nz
class SchurMat3D_right : public SchurMat3D
{
public:
    SchurMat3D_right(const int _root_nx,
                     const int _root_ny,
                     const int _root_nz,
                     const int _branch_nx,
                     const int _branch_ny,
                     const int _branch_nz)
        : SchurMat3D(_root_nx, _root_ny, _root_nz, _branch_nx, _branch_ny, _branch_nz, _root_ny * _root_nz)
    {}
    void   construct(DomainSolver3D* branch_solver) override;
    field3 operator*(const field3& root) override;
};

// Front face: y=0 plane, size = nx * nz
class SchurMat3D_front : public SchurMat3D
{
public:
    SchurMat3D_front(const int _root_nx,
                     const int _root_ny,
                     const int _root_nz,
                     const int _branch_nx,
                     const int _branch_ny,
                     const int _branch_nz)
        : SchurMat3D(_root_nx, _root_ny, _root_nz, _branch_nx, _branch_ny, _branch_nz, _root_nx * _root_nz)
    {}
    void   construct(DomainSolver3D* branch_solver) override;
    field3 operator*(const field3& root) override;
};

// Back face: y=ny-1 plane, size = nx * nz
class SchurMat3D_back : public SchurMat3D
{
public:
    SchurMat3D_back(const int _root_nx,
                    const int _root_ny,
                    const int _root_nz,
                    const int _branch_nx,
                    const int _branch_ny,
                    const int _branch_nz)
        : SchurMat3D(_root_nx, _root_ny, _root_nz, _branch_nx, _branch_ny, _branch_nz, _root_nx * _root_nz)
    {}
    void   construct(DomainSolver3D* branch_solver) override;
    field3 operator*(const field3& root) override;
};

// Down face: z=0 plane, size = nx * ny
class SchurMat3D_down : public SchurMat3D
{
public:
    SchurMat3D_down(const int _root_nx,
                    const int _root_ny,
                    const int _root_nz,
                    const int _branch_nx,
                    const int _branch_ny,
                    const int _branch_nz)
        : SchurMat3D(_root_nx, _root_ny, _root_nz, _branch_nx, _branch_ny, _branch_nz, _root_nx * _root_ny)
    {}
    void   construct(DomainSolver3D* branch_solver) override;
    field3 operator*(const field3& root) override;
};

// Up face: z=nz-1 plane, size = nx * ny
class SchurMat3D_up : public SchurMat3D
{
public:
    SchurMat3D_up(const int _root_nx,
                  const int _root_ny,
                  const int _root_nz,
                  const int _branch_nx,
                  const int _branch_ny,
                  const int _branch_nz)
        : SchurMat3D(_root_nx, _root_ny, _root_nz, _branch_nx, _branch_ny, _branch_nz, _root_nx * _root_ny)
    {}
    void   construct(DomainSolver3D* branch_solver) override;
    field3 operator*(const field3& root) override;
};
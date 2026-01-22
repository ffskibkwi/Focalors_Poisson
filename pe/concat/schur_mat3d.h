#pragma once

#include "base/pch.h"

#include "domain_solver.h"
#include "base/domain/domain3d.h"

class SchurMat3D
{
protected:
    int      cosize_n;
    int      root_nx, root_ny, root_nz, branch_nx, branch_ny, branch_nz;
    double** value;

public:
    SchurMat3D(const field3& root, const field3& branch, const int _cosize_n)
        : cosize_n(_cosize_n)
        , root_nx(root.get_nx())
        , root_ny(root.get_ny())
        , root_nz(root.get_nz())
        , branch_nx(branch.get_nx())
        , branch_ny(branch.get_ny())
        , branch_nz(branch.get_nz())
    {
        value = new double*[cosize_n];
        for (int i = 0; i < cosize_n; i++)
            value[i] = new double[cosize_n];
    }
    SchurMat3D(const Domain3DUniform& root, const Domain3DUniform& branch, const int _cosize_n)
        : cosize_n(_cosize_n)
        , root_nx(root.nx)
        , root_ny(root.ny)
        , root_nz(root.nz)
        , branch_nx(branch.nx)
        , branch_ny(branch.ny)
        , branch_nz(branch.nz)
    {
        value = new double*[cosize_n];
        for (int i = 0; i < cosize_n; i++)
            value[i] = new double[cosize_n];
    }
    ~SchurMat3D()
    {
        for (int i = 0; i < cosize_n; i++)
            delete[] value[i];
        delete[] value;
    }
    int get_size() const { return cosize_n; }
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

    virtual void   construct(DomainSolver3D* branch_solver) = 0;
    virtual field3 operator*(const field3& root)            = 0;
};

// Left face: x=0 plane, size = ny * nz
class SchurMat3D_left : public SchurMat3D
{
public:
    SchurMat3D_left(const Domain3DUniform& root, const Domain3DUniform& branch)
        : SchurMat3D(root, branch, root.ny * root.nz)
    {}
    void   construct(DomainSolver3D* branch_solver) override;
    field3 operator*(const field3& root) override;
};

// Right face: x=nx-1 plane, size = ny * nz
class SchurMat3D_right : public SchurMat3D
{
public:
    SchurMat3D_right(const Domain3DUniform& root, const Domain3DUniform& branch)
        : SchurMat3D(root, branch, root.ny * root.nz)
    {}
    void   construct(DomainSolver3D* branch_solver) override;
    field3 operator*(const field3& root) override;
};

// Front face: y=0 plane, size = nx * nz
class SchurMat3D_front : public SchurMat3D
{
public:
    SchurMat3D_front(const Domain3DUniform& root, const Domain3DUniform& branch)
        : SchurMat3D(root, branch, root.nx * root.nz)
    {}
    void   construct(DomainSolver3D* branch_solver) override;
    field3 operator*(const field3& root) override;
};

// Back face: y=ny-1 plane, size = nx * nz
class SchurMat3D_back : public SchurMat3D
{
public:
    SchurMat3D_back(const Domain3DUniform& root, const Domain3DUniform& branch)
        : SchurMat3D(root, branch, root.nx * root.nz)
    {}
    void   construct(DomainSolver3D* branch_solver) override;
    field3 operator*(const field3& root) override;
};

// Down face: z=0 plane, size = nx * ny
class SchurMat3D_down : public SchurMat3D
{
public:
    SchurMat3D_down(const Domain3DUniform& root, const Domain3DUniform& branch)
        : SchurMat3D(root, branch, root.nx * root.ny)
    {}
    void   construct(DomainSolver3D* branch_solver) override;
    field3 operator*(const field3& root) override;
};

// Up face: z=nz-1 plane, size = nx * ny
class SchurMat3D_up : public SchurMat3D
{
public:
    SchurMat3D_up(const Domain3DUniform& root, const Domain3DUniform& branch)
        : SchurMat3D(root, branch, root.nx * root.ny)
    {}
    void   construct(DomainSolver3D* branch_solver) override;
    field3 operator*(const field3& root) override;
};
#pragma once

#include "base/pch.h"

#include "domain_solver.h"
#include "base/domain/domain3d.h"

class Schur_mat3d
{
protected:
    int      cosize_n;
    int      root_nx, root_ny, root_nz, branch_nx, branch_ny, branch_nz;
    double** value;

public:
    Schur_mat3d(const field3& root, const field3& branch, const int _cosize_n)
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
    Schur_mat3d(const Domain3DUniform& root, const Domain3DUniform& branch, const int _cosize_n)
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
    ~Schur_mat3d()
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
class Schur_mat3d_left : public Schur_mat3d
{
public:
    Schur_mat3d_left(const Domain3DUniform& root, const Domain3DUniform& branch)
        : Schur_mat3d(root, branch, root.ny * root.nz)
    {}
    void   construct(DomainSolver3D* branch_solver) override;
    field3 operator*(const field3& root) override;
};

// Right face: x=nx-1 plane, size = ny * nz
class Schur_mat3d_right : public Schur_mat3d
{
public:
    Schur_mat3d_right(const Domain3DUniform& root, const Domain3DUniform& branch)
        : Schur_mat3d(root, branch, root.ny * root.nz)
    {}
    void   construct(DomainSolver3D* branch_solver) override;
    field3 operator*(const field3& root) override;
};

// Front face: y=0 plane, size = nx * nz
class Schur_mat3d_front : public Schur_mat3d
{
public:
    Schur_mat3d_front(const Domain3DUniform& root, const Domain3DUniform& branch)
        : Schur_mat3d(root, branch, root.nx * root.nz)
    {}
    void   construct(DomainSolver3D* branch_solver) override;
    field3 operator*(const field3& root) override;
};

// Back face: y=ny-1 plane, size = nx * nz
class Schur_mat3d_back : public Schur_mat3d
{
public:
    Schur_mat3d_back(const Domain3DUniform& root, const Domain3DUniform& branch)
        : Schur_mat3d(root, branch, root.nx * root.nz)
    {}
    void   construct(DomainSolver3D* branch_solver) override;
    field3 operator*(const field3& root) override;
};

// Down face: z=0 plane, size = nx * ny
class Schur_mat3d_down : public Schur_mat3d
{
public:
    Schur_mat3d_down(const Domain3DUniform& root, const Domain3DUniform& branch)
        : Schur_mat3d(root, branch, root.nx * root.ny)
    {}
    void   construct(DomainSolver3D* branch_solver) override;
    field3 operator*(const field3& root) override;
};

// Up face: z=nz-1 plane, size = nx * ny
class Schur_mat3d_up : public Schur_mat3d
{
public:
    Schur_mat3d_up(const Domain3DUniform& root, const Domain3DUniform& branch)
        : Schur_mat3d(root, branch, root.nx * root.ny)
    {}
    void   construct(DomainSolver3D* branch_solver) override;
    field3 operator*(const field3& root) override;
};
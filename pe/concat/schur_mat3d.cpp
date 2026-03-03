#include "schur_mat3d.h"

#include "domain_solver.h"

// XNegative: 处理来自左侧邻域的贡献（作用于自身的 x=0 界面）
void SchurMat3D_xneg::construct(DomainSolver3D* branch_solver)
{
    int    interface_size = bny * bnz;
    field3 t_a(bnx, bny, bnz);
    for (int j = 0; j < bny; j++)
    {
        for (int k = 0; k < bnz; k++)
        {
            int col_idx = j * bnz + k;
            t_a.clear();
            t_a(bnx - 1, j, k) = 1.0;
            branch_solver->solve(t_a);
            for (int jj = 0; jj < bny; jj++)
            {
                for (int kk = 0; kk < bnz; kk++)
                {
                    int row_idx             = jj * bnz + kk;
                    value(row_idx, col_idx) = t_a(bnx - 1, jj, kk);
                }
            }
        }
    }
}

field3 SchurMat3D_xneg::operator*(const field3& root)
{
    int root_nx = root.get_nx();
    int root_ny = root.get_ny();
    int root_nz = root.get_nz();

    field3 R(root_nx, root_ny, root_nz);
    OPENMP_PARALLEL_FOR()
    for (int i = 0; i < root_ny * root_nz; i++)
    {
        int iy       = i / root_nz;
        int iz       = i % root_nz;
        R(0, iy, iz) = 0.;
        for (int j = 0; j < root_ny * root_nz; j++)
        {
            int jy = j / root_nz;
            int jz = j % root_nz;
            R(0, iy, iz) += root(0, jy, jz) * value(i, j);
        }
    }
    return R;
}

// XPositive: 处理来自右侧邻域的贡献（作用于自身的 x=nx-1 界面）
void SchurMat3D_xpos::construct(DomainSolver3D* branch_solver)
{
    field3 t_a(bnx, bny, bnz);
    for (int j = 0; j < bny; j++)
    {
        for (int k = 0; k < bnz; k++)
        {
            int col_idx = j * bnz + k;
            t_a.clear();
            t_a(0, j, k) = 1.0;
            branch_solver->solve(t_a);
            for (int jj = 0; jj < bny; jj++)
            {
                for (int kk = 0; kk < bnz; kk++)
                {
                    int row_idx             = jj * bnz + kk;
                    value(row_idx, col_idx) = t_a(0, jj, kk);
                }
            }
        }
    }
}

field3 SchurMat3D_xpos::operator*(const field3& root)
{
    int root_nx = root.get_nx();
    int root_ny = root.get_ny();
    int root_nz = root.get_nz();

    field3 R(root_nx, root_ny, root_nz);
    OPENMP_PARALLEL_FOR()
    for (int i = 0; i < root_ny * root_nz; i++)
    {
        int iy                 = i / root_nz;
        int iz                 = i % root_nz;
        R(root_nx - 1, iy, iz) = 0.;
        for (int j = 0; j < root_ny * root_nz; j++)
        {
            int jy = j / root_nz;
            int jz = j % root_nz;
            R(root_nx - 1, iy, iz) += root(root_nx - 1, jy, jz) * value(i, j);
        }
    }
    return R;
}

// YNegative (y 负方向): 作用于自身的 y=0 界面
void SchurMat3D_yneg::construct(DomainSolver3D* branch_solver)
{
    field3 t_a(bnx, bny, bnz);
    for (int i = 0; i < bnx; i++)
    {
        for (int k = 0; k < bnz; k++)
        {
            int col_idx = i * bnz + k;
            t_a.clear();
            t_a(i, bny - 1, k) = 1.0;
            branch_solver->solve(t_a);
            for (int ii = 0; ii < bnx; ii++)
            {
                for (int kk = 0; kk < bnz; kk++)
                {
                    int row_idx             = ii * bnz + kk;
                    value(row_idx, col_idx) = t_a(ii, bny - 1, kk);
                }
            }
        }
    }
}

field3 SchurMat3D_yneg::operator*(const field3& root)
{
    int root_nx = root.get_nx();
    int root_ny = root.get_ny();
    int root_nz = root.get_nz();

    field3 R(root_nx, root_ny, root_nz);
    OPENMP_PARALLEL_FOR()
    for (int idx = 0; idx < root_nx * root_nz; idx++)
    {
        int ix       = idx / root_nz;
        int iz       = idx % root_nz;
        R(ix, 0, iz) = 0.;
        for (int jdx = 0; jdx < root_nx * root_nz; jdx++)
        {
            int jx = jdx / root_nz;
            int jz = jdx % root_nz;
            R(ix, 0, iz) += root(jx, 0, jz) * value(idx, jdx);
        }
    }
    return R;
}

// YPositive (y 正方向): 作用于自身的 y=ny-1 界面
void SchurMat3D_ypos::construct(DomainSolver3D* branch_solver)
{
    field3 t_a(bnx, bny, bnz);
    for (int i = 0; i < bnx; i++)
    {
        for (int k = 0; k < bnz; k++)
        {
            int col_idx = i * bnz + k;
            t_a.clear();
            t_a(i, 0, k) = 1.0;
            branch_solver->solve(t_a);
            for (int ii = 0; ii < bnx; ii++)
            {
                for (int kk = 0; kk < bnz; kk++)
                {
                    int row_idx             = ii * bnz + kk;
                    value(row_idx, col_idx) = t_a(ii, 0, kk);
                }
            }
        }
    }
}

field3 SchurMat3D_ypos::operator*(const field3& root)
{
    int root_nx = root.get_nx();
    int root_ny = root.get_ny();
    int root_nz = root.get_nz();

    field3 R(root_nx, root_ny, root_nz);
    OPENMP_PARALLEL_FOR()
    for (int idx = 0; idx < root_nx * root_nz; idx++)
    {
        int ix                 = idx / root_nz;
        int iz                 = idx % root_nz;
        R(ix, root_ny - 1, iz) = 0.;
        for (int jdx = 0; jdx < root_nx * root_nz; jdx++)
        {
            int jx = jdx / root_nz;
            int jz = jdx % root_nz;
            R(ix, root_ny - 1, iz) += root(jx, root_ny - 1, jz) * value(idx, jdx);
        }
    }
    return R;
}

// ZNegative (z 负方向): 作用于自身的 z=0 界面
void SchurMat3D_zneg::construct(DomainSolver3D* branch_solver)
{
    field3 t_a(bnx, bny, bnz);
    for (int i = 0; i < bnx; i++)
    {
        for (int j = 0; j < bny; j++)
        {
            int col_idx = i * bny + j;
            t_a.clear();
            t_a(i, j, bnz - 1) = 1.0;
            branch_solver->solve(t_a);
            for (int ii = 0; ii < bnx; ii++)
            {
                for (int jj = 0; jj < bny; jj++)
                {
                    int row_idx             = ii * bny + jj;
                    value(row_idx, col_idx) = t_a(ii, jj, bnz - 1);
                }
            }
        }
    }
}

field3 SchurMat3D_zneg::operator*(const field3& root)
{
    int root_nx = root.get_nx();
    int root_ny = root.get_ny();
    int root_nz = root.get_nz();

    field3 R(root_nx, root_ny, root_nz);
    OPENMP_PARALLEL_FOR()
    for (int idx = 0; idx < root_nx * root_ny; idx++)
    {
        int ix       = idx / root_ny;
        int iy       = idx % root_ny;
        R(ix, iy, 0) = 0.;
        for (int jdx = 0; jdx < root_nx * root_ny; jdx++)
        {
            int jx = jdx / root_ny;
            int jy = jdx % root_ny;
            R(ix, iy, 0) += root(jx, jy, 0) * value(idx, jdx);
        }
    }
    return R;
}

// ZPositive (z 正方向): 作用于自身的 z=nz-1 界面
void SchurMat3D_zpos::construct(DomainSolver3D* branch_solver)
{
    field3 t_a(bnx, bny, bnz);
    for (int i = 0; i < bnx; i++)
    {
        for (int j = 0; j < bny; j++)
        {
            int col_idx = i * bny + j;
            t_a.clear();
            t_a(i, j, 0) = 1.0;
            branch_solver->solve(t_a);
            for (int ii = 0; ii < bnx; ii++)
            {
                for (int jj = 0; jj < bny; jj++)
                {
                    int row_idx             = ii * bny + jj;
                    value(row_idx, col_idx) = t_a(ii, jj, 0);
                }
            }
        }
    }
}

field3 SchurMat3D_zpos::operator*(const field3& root)
{
    int root_nx = root.get_nx();
    int root_ny = root.get_ny();
    int root_nz = root.get_nz();

    field3 R(root_nx, root_ny, root_nz);
    OPENMP_PARALLEL_FOR()
    for (int idx = 0; idx < root_nx * root_ny; idx++)
    {
        int ix                 = idx / root_ny;
        int iy                 = idx % root_ny;
        R(ix, iy, root_nz - 1) = 0.;
        for (int jdx = 0; jdx < root_nx * root_ny; jdx++)
        {
            int jx = jdx / root_ny;
            int jy = jdx % root_ny;
            R(ix, iy, root_nz - 1) += root(jx, jy, root_nz - 1) * value(idx, jdx);
        }
    }
    return R;
}
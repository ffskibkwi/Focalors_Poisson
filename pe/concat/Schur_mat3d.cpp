#include "Schur_mat3d.h"

void Schur_mat3d_left::construct(DomainSolver3D* branch_solver)
{
    field3 t_a(branch_nx, branch_ny, branch_nz, "temp");
    int idx = 0;
    for (int j = 0; j < root_ny; j++)
    {
        for (int k = 0; k < root_nz; k++)
        {
            t_a.clear();
            t_a(branch_nx - 1, j, k) = 1.;
            branch_solver->solve(t_a);
            int col_idx = 0;
            for (int j2 = 0; j2 < root_ny; j2++)
            {
                for (int k2 = 0; k2 < root_nz; k2++)
                {
                    value[col_idx][idx] = t_a(branch_nx - 1, j2, k2);
                    col_idx++;
                }
            }
            idx++;
        }
    }
}

field3 Schur_mat3d_left::operator*(const field3& root)
{
    field3 R(root_nx, root_ny, root_nz, "result");

    OPENMP_PARALLEL_FOR()
    for (int j = 0; j < root_ny; j++)
    {
        for (int k = 0; k < root_nz; k++)
        {
            int row_idx = j * root_nz + k;
            R(0, j, k) = 0.;
            int col_idx = 0;
            for (int j2 = 0; j2 < root_ny; j2++)
            {
                for (int k2 = 0; k2 < root_nz; k2++)
                {
                    R(0, j, k) += root(0, j2, k2) * value[row_idx][col_idx];
                    col_idx++;
                }
            }
        }
    }
    return R;
}

void Schur_mat3d_right::construct(DomainSolver3D* branch_solver)
{
    field3 t_a(branch_nx, branch_ny, branch_nz, "temp");
    int idx = 0;
    for (int j = 0; j < root_ny; j++)
    {
        for (int k = 0; k < root_nz; k++)
        {
            t_a.clear();
            t_a(0, j, k) = 1.;
            branch_solver->solve(t_a);
            int col_idx = 0;
            for (int j2 = 0; j2 < root_ny; j2++)
            {
                for (int k2 = 0; k2 < root_nz; k2++)
                {
                    value[col_idx][idx] = t_a(0, j2, k2);
                    col_idx++;
                }
            }
            idx++;
        }
    }
}

field3 Schur_mat3d_right::operator*(const field3& root)
{
    field3 R(root_nx, root_ny, root_nz, "result");

    OPENMP_PARALLEL_FOR()
    for (int j = 0; j < root_ny; j++)
    {
        for (int k = 0; k < root_nz; k++)
        {
            int row_idx = j * root_nz + k;
            R(root_nx - 1, j, k) = 0.;
            int col_idx = 0;
            for (int j2 = 0; j2 < root_ny; j2++)
            {
                for (int k2 = 0; k2 < root_nz; k2++)
                {
                    R(root_nx - 1, j, k) += root(root_nx - 1, j2, k2) * value[row_idx][col_idx];
                    col_idx++;
                }
            }
        }
    }
    return R;
}

void Schur_mat3d_front::construct(DomainSolver3D* branch_solver)
{
    field3 t_a(branch_nx, branch_ny, branch_nz, "temp");
    int idx = 0;
    for (int i = 0; i < root_nx; i++)
    {
        for (int k = 0; k < root_nz; k++)
        {
            t_a.clear();
            t_a(i, branch_ny - 1, k) = 1.;
            branch_solver->solve(t_a);
            int col_idx = 0;
            for (int i2 = 0; i2 < root_nx; i2++)
            {
                for (int k2 = 0; k2 < root_nz; k2++)
                {
                    value[col_idx][idx] = t_a(i2, branch_ny - 1, k2);
                    col_idx++;
                }
            }
            idx++;
        }
    }
}

field3 Schur_mat3d_front::operator*(const field3& root)
{
    field3 R(root_nx, root_ny, root_nz, "result");

    OPENMP_PARALLEL_FOR()
    for (int i = 0; i < root_nx; i++)
    {
        for (int k = 0; k < root_nz; k++)
        {
            int row_idx = i * root_nz + k;
            R(i, 0, k) = 0.;
            int col_idx = 0;
            for (int i2 = 0; i2 < root_nx; i2++)
            {
                for (int k2 = 0; k2 < root_nz; k2++)
                {
                    R(i, 0, k) += root(i2, 0, k2) * value[row_idx][col_idx];
                    col_idx++;
                }
            }
        }
    }
    return R;
}

void Schur_mat3d_back::construct(DomainSolver3D* branch_solver)
{
    field3 t_a(branch_nx, branch_ny, branch_nz, "temp");
    int idx = 0;
    for (int i = 0; i < root_nx; i++)
    {
        for (int k = 0; k < root_nz; k++)
        {
            t_a.clear();
            t_a(i, 0, k) = 1.;
            branch_solver->solve(t_a);
            int col_idx = 0;
            for (int i2 = 0; i2 < root_nx; i2++)
            {
                for (int k2 = 0; k2 < root_nz; k2++)
                {
                    value[col_idx][idx] = t_a(i2, 0, k2);
                    col_idx++;
                }
            }
            idx++;
        }
    }
}

field3 Schur_mat3d_back::operator*(const field3& root)
{
    field3 R(root_nx, root_ny, root_nz, "result");

    OPENMP_PARALLEL_FOR()
    for (int i = 0; i < root_nx; i++)
    {
        for (int k = 0; k < root_nz; k++)
        {
            int row_idx = i * root_nz + k;
            R(i, root_ny - 1, k) = 0.;
            int col_idx = 0;
            for (int i2 = 0; i2 < root_nx; i2++)
            {
                for (int k2 = 0; k2 < root_nz; k2++)
                {
                    R(i, root_ny - 1, k) += root(i2, root_ny - 1, k2) * value[row_idx][col_idx];
                    col_idx++;
                }
            }
        }
    }
    return R;
}

void Schur_mat3d_down::construct(DomainSolver3D* branch_solver)
{
    field3 t_a(branch_nx, branch_ny, branch_nz, "temp");
    int idx = 0;
    for (int i = 0; i < root_nx; i++)
    {
        for (int j = 0; j < root_ny; j++)
        {
            t_a.clear();
            t_a(i, j, branch_nz - 1) = 1.;
            branch_solver->solve(t_a);
            int col_idx = 0;
            for (int i2 = 0; i2 < root_nx; i2++)
            {
                for (int j2 = 0; j2 < root_ny; j2++)
                {
                    value[col_idx][idx] = t_a(i2, j2, branch_nz - 1);
                    col_idx++;
                }
            }
            idx++;
        }
    }
}

field3 Schur_mat3d_down::operator*(const field3& root)
{
    field3 R(root_nx, root_ny, root_nz, "result");

    OPENMP_PARALLEL_FOR()
    for (int i = 0; i < root_nx; i++)
    {
        for (int j = 0; j < root_ny; j++)
        {
            int row_idx = i * root_ny + j;
            R(i, j, 0) = 0.;
            int col_idx = 0;
            for (int i2 = 0; i2 < root_nx; i2++)
            {
                for (int j2 = 0; j2 < root_ny; j2++)
                {
                    R(i, j, 0) += root(i2, j2, 0) * value[row_idx][col_idx];
                    col_idx++;
                }
            }
        }
    }
    return R;
}

void Schur_mat3d_up::construct(DomainSolver3D* branch_solver)
{
    field3 t_a(branch_nx, branch_ny, branch_nz, "temp");
    int idx = 0;
    for (int i = 0; i < root_nx; i++)
    {
        for (int j = 0; j < root_ny; j++)
        {
            t_a.clear();
            t_a(i, j, 0) = 1.;
            branch_solver->solve(t_a);
            int col_idx = 0;
            for (int i2 = 0; i2 < root_nx; i2++)
            {
                for (int j2 = 0; j2 < root_ny; j2++)
                {
                    value[col_idx][idx] = t_a(i2, j2, 0);
                    col_idx++;
                }
            }
            idx++;
        }
    }
}

field3 Schur_mat3d_up::operator*(const field3& root)
{
    field3 R(root_nx, root_ny, root_nz, "result");

    OPENMP_PARALLEL_FOR()
    for (int i = 0; i < root_nx; i++)
    {
        for (int j = 0; j < root_ny; j++)
        {
            int row_idx = i * root_ny + j;
            R(i, j, root_nz - 1) = 0.;
            int col_idx = 0;
            for (int i2 = 0; i2 < root_nx; i2++)
            {
                for (int j2 = 0; j2 < root_ny; j2++)
                {
                    R(i, j, root_nz - 1) += root(i2, j2, root_nz - 1) * value[row_idx][col_idx];
                    col_idx++;
                }
            }
        }
    }
    return R;
}
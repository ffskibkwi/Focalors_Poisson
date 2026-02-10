#include "schur_mat2d.h"

#include "domain_solver.h"

void SchurMat2D_left::construct(DomainSolver2D* branch_solver)
{
    field2 t_a(bnx, bny);
    for (int i = 0; i < cn; i++)
    {
        t_a.clear();
        t_a(bnx - 1, i) = 1.;
        branch_solver->solve(t_a);
        for (int j = 0; j < cn; j++)
            value(j, i) = t_a(bnx - 1, j);
    }
}

field2 SchurMat2D_left::operator*(const field2& root)
{
    int rnx = root.get_nx();
    int rny = root.get_ny();

    field2 R(rnx, rny);

    OPENMP_PARALLEL_FOR()
    for (int i = 0; i < cn; i++)
    {
        R(0, i) = 0.;
        for (int j = 0; j < cn; j++)
            R(0, i) += root(0, j) * value(i, j);
    }
    return R;
}

void SchurMat2D_right::construct(DomainSolver2D* branch_solver)
{
    field2 t_a(bnx, bny);
    for (int i = 0; i < cn; i++)
    {
        t_a.clear();
        t_a(0, i) = 1.;
        branch_solver->solve(t_a);
        for (int j = 0; j < cn; j++)
            value(j, i) = t_a(0, j);
    }
}

field2 SchurMat2D_right::operator*(const field2& root)
{
    int rnx = root.get_nx();
    int rny = root.get_ny();

    field2 R(rnx, rny);

    OPENMP_PARALLEL_FOR()
    for (int i = 0; i < cn; i++)
    {
        R(rnx - 1, i) = 0.;
        for (int j = 0; j < cn; j++)
            R(rnx - 1, i) += root(rnx - 1, j) * value(i, j);
    }
    return R;
}

void SchurMat2D_up::construct(DomainSolver2D* branch_solver)
{
    field2 t_a(bnx, bny);
    for (int i = 0; i < cn; i++)
    {
        t_a.clear();
        t_a(i, 0) = 1.;
        branch_solver->solve(t_a);
        for (int j = 0; j < cn; j++)
            value(j, i) = t_a(j, 0);
    }
}

field2 SchurMat2D_up::operator*(const field2& root)
{
    int rnx = root.get_nx();
    int rny = root.get_ny();

    field2 R(rnx, rny);

    OPENMP_PARALLEL_FOR()
    for (int i = 0; i < cn; i++)
    {
        R(i, rny - 1) = 0.;
        for (int j = 0; j < cn; j++)
            R(i, rny - 1) += root(j, rny - 1) * value(i, j);
    }
    return R;
}

void SchurMat2D_down::construct(DomainSolver2D* branch_solver)
{
    field2 t_a(bnx, bny);
    for (int i = 0; i < cn; i++)
    {
        t_a.clear();
        t_a(i, bny - 1) = 1.;
        branch_solver->solve(t_a);
        for (int j = 0; j < cn; j++)
            value(j, i) = t_a(j, bny - 1);
    }
}

field2 SchurMat2D_down::operator*(const field2& root)
{
    int rnx = root.get_nx();
    int rny = root.get_ny();

    field2 R(rnx, rny);

    OPENMP_PARALLEL_FOR()
    for (int i = 0; i < cn; i++)
    {
        R(i, 0) = 0.;
        for (int j = 0; j < cn; j++)
            R(i, 0) += root(j, 0) * value(i, j);
    }
    return R;
}
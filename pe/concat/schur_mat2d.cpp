#include "schur_mat2d.h"

#include "domain_solver.h"
#include "io/csv_writer_2d.h"

void SchurMat2D::dump_to_csv(const std::string& directory)
{
    std::string filename = directory + "/" + name + ".csv";
    // Check if IO module available
    IO::array_to_csv(value, concat_n, concat_n, filename);
}

void SchurMat2D_left::construct(DomainSolver2D* branch_solver)
{
    field2 t_a(branch_nx, branch_ny);
    for (int i = 0; i < concat_n; i++)
    {
        t_a.clear();
        t_a(branch_nx - 1, i) = 1.;
        branch_solver->solve(t_a, false);
        for (int j = 0; j < concat_n; j++)
            value[j][i] = t_a(branch_nx - 1, j);
    }
}

field2 SchurMat2D_left::operator*(const field2& root)
{
    field2 R(root.get_nx(), root.get_ny());

    OPENMP_PARALLEL_FOR()
    for (int i = 0; i < concat_n; i++)
    {
        R(0, i) = 0.;
        for (int j = 0; j < concat_n; j++)
            R(0, i) += root(0, j) * value[i][j];
    }
    return R;
}

void SchurMat2D_right::construct(DomainSolver2D* branch_solver)
{
    field2 t_a(branch_nx, branch_ny);
    for (int i = 0; i < concat_n; i++)
    {
        t_a.clear();
        t_a(0, i) = 1.;
        branch_solver->solve(t_a, false);
        for (int j = 0; j < concat_n; j++)
            value[j][i] = t_a(0, j);
    }
}

field2 SchurMat2D_right::operator*(const field2& root)
{
    field2 R(root.get_nx(), root.get_ny());

    int root_nx = root.get_nx();

    OPENMP_PARALLEL_FOR()
    for (int i = 0; i < concat_n; i++)
    {
        R(root_nx - 1, i) = 0.;
        for (int j = 0; j < concat_n; j++)
            R(root_nx - 1, i) += root(root_nx - 1, j) * value[i][j];
    }
    return R;
}

void SchurMat2D_up::construct(DomainSolver2D* branch_solver)
{
    field2 t_a(branch_nx, branch_ny);
    for (int i = 0; i < concat_n; i++)
    {
        t_a.clear();
        t_a(i, 0) = 1.;
        branch_solver->solve(t_a, false);
        for (int j = 0; j < concat_n; j++)
            value[j][i] = t_a(j, 0);
    }
}

field2 SchurMat2D_up::operator*(const field2& root)
{
    field2 R(root.get_nx(), root.get_ny());

    int root_ny = root.get_nx();

    OPENMP_PARALLEL_FOR()
    for (int i = 0; i < concat_n; i++)
    {
        R(i, root_ny - 1) = 0.;
        for (int j = 0; j < concat_n; j++)
            R(i, root_ny - 1) += root(j, root_ny - 1) * value[i][j];
    }
    return R;
}

void SchurMat2D_down::construct(DomainSolver2D* branch_solver)
{
    field2 t_a(branch_nx, branch_ny);
    for (int i = 0; i < concat_n; i++)
    {
        t_a.clear();
        t_a(i, branch_ny - 1) = 1.;
        branch_solver->solve(t_a, false);
        for (int j = 0; j < concat_n; j++)
            value[j][i] = t_a(j, branch_ny - 1);
    }
}

field2 SchurMat2D_down::operator*(const field2& root)
{
    field2 R(root.get_nx(), root.get_ny());

    OPENMP_PARALLEL_FOR()
    for (int i = 0; i < concat_n; i++)
    {
        R(i, 0) = 0.;
        for (int j = 0; j < concat_n; j++)
            R(i, 0) += root(j, 0) * value[i][j];
    }
    return R;
}
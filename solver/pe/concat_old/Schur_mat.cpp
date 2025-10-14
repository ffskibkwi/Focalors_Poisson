#include "Schur_mat.h"

Schur_mat::Schur_mat(const field2& root, const field2& branch, LocationType dir)
    : direction(dir)
{
    switch (dir)
    {
    case LocationType::Left:
    case LocationType::Right:
        cosize_n = root.get_ny();
        break;
    case LocationType::Up:
    case LocationType::Down:
        cosize_n = root.get_nx();
        break;
    default:
        throw std::invalid_argument("Schur_mat only supports 2D directions");
    }

    root_nx   = root.get_nx();
    root_ny   = root.get_ny();
    branch_nx = branch.get_nx();
    branch_ny = branch.get_ny();

    value = new double*[cosize_n];
    for (int i = 0; i < cosize_n; i++)
        value[i] = new double[cosize_n];
}

Schur_mat::~Schur_mat()
{
    for (int i = 0; i < cosize_n; i++)
        delete[] value[i];
    delete[] value;
}

void Schur_mat::print()
{
    for (int i = 0; i < cosize_n; i++)
    {
        for (int j = 0; j < cosize_n; j++)
            std::cout << value[i][j] << std::endl;
    }
}

void Schur_mat::construct(PoissonSolver2DInterface& branch_solver)
{
    field2 t_a(branch_nx, branch_ny);

    switch (direction)
    {
    case LocationType::Left:
        for (int i = 0; i < cosize_n; i++)
        {
            t_a.clear();
            t_a(branch_nx - 1, i) = 1.;
            branch_solver.solve(t_a);
            for (int j = 0; j < cosize_n; j++)
                value[j][i] = t_a(branch_nx - 1, j);
        }
        break;
    case LocationType::Right:
        for (int i = 0; i < cosize_n; i++)
        {
            t_a.clear();
            t_a(0, i) = 1.;
            branch_solver.solve(t_a);
            for (int j = 0; j < cosize_n; j++)
                value[j][i] = t_a(0, j);
        }
        break;
    case LocationType::Up:
        for (int i = 0; i < cosize_n; i++)
        {
            t_a.clear();
            t_a(i, 0) = 1.;
            branch_solver.solve(t_a);
            for (int j = 0; j < cosize_n; j++)
                value[j][i] = t_a(j, 0);
        }
        break;
    case LocationType::Down:
        for (int i = 0; i < cosize_n; i++)
        {
            t_a.clear();
            t_a(i, branch_ny - 1) = 1.;
            branch_solver.solve(t_a);
            for (int j = 0; j < cosize_n; j++)
                value[j][i] = t_a(j, branch_ny - 1);
        }
        break;
    default:
        throw std::invalid_argument("Schur_mat only supports 2D directions");
    }
}

field2 Schur_mat::operator*(const field2& root)
{
    field2 R(root_nx, root_ny);

    switch (direction)
    {
    case LocationType::Left:
        OPENMP_PARALLEL_FOR()
        for (int i = 0; i < cosize_n; i++)
        {
            R(0, i) = 0.;
            for (int j = 0; j < cosize_n; j++)
                R(0, i) += root(0, j) * value[i][j];
        }
        break;
    case LocationType::Right:
        OPENMP_PARALLEL_FOR()
        for (int i = 0; i < cosize_n; i++)
        {
            R(root_nx - 1, i) = 0.;
            for (int j = 0; j < cosize_n; j++)
                R(root_nx - 1, i) += root(root_nx - 1, j) * value[i][j];
        }
        break;
    case LocationType::Up:
        OPENMP_PARALLEL_FOR()
        for (int i = 0; i < cosize_n; i++)
        {
            R(i, root_ny - 1) = 0.;
            for (int j = 0; j < cosize_n; j++)
                R(i, root_ny - 1) += root(i, root_ny - 1) * value[i][j];
        }
        break;
    case LocationType::Down:
        OPENMP_PARALLEL_FOR()
        for (int i = 0; i < cosize_n; i++)
        {
            R(i, 0) = 0.;
            for (int j = 0; j < cosize_n; j++)
                R(i, 0) += root(i, 0) * value[i][j];
        }
        break;
    default:
        throw std::invalid_argument("Schur_mat only supports 2D directions");
    }

    return R;
}



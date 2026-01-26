#include "schur_mat2d_slab_x.h"

#include "domain_solver.h"
#include "io/csv_writer_2d.h"

void SchurMat2DSlabX_left::construct(DomainSolver2D* branch_solver)
{
    field2 t_a(bsnx, bny);
    for (int i = 0; i < cn; i++)
    {
        t_a.clear();
        if (mpi_rank == mpi_size - 1)
            t_a(bsnx - 1, i) = 1.;
        branch_solver->solve(t_a, false);
        MPI_Scatterv(t_a.get_ptr(bsnx - 1),
                     cs_lengths,
                     cs_displacements,
                     MPI_DOUBLE,
                     buf_csn,
                     csn,
                     MPI_DOUBLE,
                     mpi_size - 1,
                     communicator);
        for (int j = 0; j < csn; j++)
            value(j, i) = buf_csn[j];
    }
}

field2 SchurMat2DSlabX_left::operator*(const field2& root)
{
    int rsnx = root.get_nx(); // root domain with slab decomposition nx
    int rny  = root.get_ny();

    field2 R(rsnx, rny);

    if (mpi_rank == 0)
    {
        for (int j = 0; j < cn; j++)
            buf_cn[j] = root(0, j);
    }
    MPI_Bcast(buf_cn, cn, MPI_DOUBLE, 0, communicator);

    OPENMP_PARALLEL_FOR()
    for (int i = 0; i < csn; i++)
    {
        buf_csn[i] = 0.;
        for (int j = 0; j < cn; j++)
            buf_csn[i] += buf_cn[j] * value(i, j);
    }
    MPI_Gatherv(buf_csn, csn, MPI_DOUBLE, R.get_ptr(0), cs_lengths, cs_displacements, MPI_DOUBLE, 0, communicator);
    return R;
}

void SchurMat2DSlabX_right::construct(DomainSolver2D* branch_solver)
{
    field2 t_a(bnx, bny);
    for (int i = 0; i < cn; i++)
    {
        t_a.clear();
        t_a(0, i) = 1.;
        branch_solver->solve(t_a, false);
        for (int j = 0; j < cn; j++)
            value(j, i) = t_a(0, j);
    }
}

field2 SchurMat2DSlabX_right::operator*(const field2& root)
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

void SchurMat2DSlabX_up::construct(DomainSolver2D* branch_solver)
{
    field2 t_a(bnx, bny);
    for (int i = 0; i < cn; i++)
    {
        t_a.clear();
        t_a(i, 0) = 1.;
        branch_solver->solve(t_a, false);
        for (int j = 0; j < cn; j++)
            value(j, i) = t_a(j, 0);
    }
}

field2 SchurMat2DSlabX_up::operator*(const field2& root)
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

void SchurMat2DSlabX_down::construct(DomainSolver2D* branch_solver)
{
    field2 t_a(bnx, bny);
    for (int i = 0; i < cn; i++)
    {
        t_a.clear();
        t_a(i, bny - 1) = 1.;
        branch_solver->solve(t_a, false);
        for (int j = 0; j < cn; j++)
            value(j, i) = t_a(j, bny - 1);
    }
}

field2 SchurMat2DSlabX_down::operator*(const field2& root)
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
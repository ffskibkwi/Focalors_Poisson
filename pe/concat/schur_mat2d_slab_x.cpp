#include "schur_mat2d_slab_x.h"

#include "base/parallel/mpi/distribute_slab.h"
#include "domain_solver.h"
#include "io/csv_writer_2d.h"

void SchurMat2DSlabX::write_csv(const std::string& directory) { IO::write_csv(value, directory); }

void migrate_from(SchurMat2DSlabX* src, SchurMat2DSlabX* dest)
{
    field2*  data_src  = src != nullptr ? &(src->value) : nullptr;
    field2*  data_dest = dest != nullptr ? &(dest->value) : nullptr;
    MPI_Comm comm_src  = src != nullptr ? src->communicator : MPI_COMM_NULL;
    MPI_Comm comm_dest = dest != nullptr ? dest->communicator : MPI_COMM_NULL;
    MPIUtils::redistribute_slab(data_src, data_dest, comm_src, comm_dest);
}

void SchurMat2DSlabX_left::construct(DomainSolver2D* branch_solver)
{
    field2 t_a(bsnx, bny);
    for (int i = 0; i < cn; i++)
    {
        t_a.clear();
        if (mpi_rank == mpi_size - 1)
            t_a(bsnx - 1, i) = 1.;
        branch_solver->solve(t_a);
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
    field2 t_a(bsnx, bny);
    for (int i = 0; i < cn; i++)
    {
        t_a.clear();
        if (mpi_rank == 0)
            t_a(0, i) = 1.;
        branch_solver->solve(t_a);
        MPI_Scatterv(
            t_a.get_ptr(0), cs_lengths, cs_displacements, MPI_DOUBLE, buf_csn, csn, MPI_DOUBLE, 0, communicator);
        for (int j = 0; j < csn; j++)
            value(j, i) = buf_csn[j];
    }
}

field2 SchurMat2DSlabX_right::operator*(const field2& root)
{
    int rsnx = root.get_nx();
    int rny  = root.get_ny();

    field2 R(rsnx, rny);

    if (mpi_rank == mpi_size - 1)
    {
        for (int j = 0; j < cn; j++)
            buf_cn[j] = root(rsnx - 1, j);
    }
    MPI_Bcast(buf_cn, cn, MPI_DOUBLE, mpi_size - 1, communicator);

    OPENMP_PARALLEL_FOR()
    for (int i = 0; i < csn; i++)
    {
        buf_csn[i] = 0.;
        for (int j = 0; j < cn; j++)
            buf_csn[i] += buf_cn[j] * value(i, j);
    }
    MPI_Gatherv(buf_csn,
                csn,
                MPI_DOUBLE,
                R.get_ptr(rsnx - 1),
                cs_lengths,
                cs_displacements,
                MPI_DOUBLE,
                mpi_size - 1,
                communicator);
    return R;
}

void SchurMat2DSlabX_up::construct(DomainSolver2D* branch_solver)
{
    field2 t_a(bsnx, bny);
    for (int i = 0; i < cn; i++)
    {
        t_a.clear();
        if (i >= cs_displacements[mpi_rank])
        {
            if (mpi_rank < mpi_size - 1 && i >= cs_displacements[mpi_rank + 1])
                continue;
            t_a(i - cs_displacements[mpi_rank], 0) = 1.;
        }
        branch_solver->solve(t_a);
        for (int j = 0; j < csn; j++)
            value(j, i) = t_a(j, 0);
    }
}

field2 SchurMat2DSlabX_up::operator*(const field2& root)
{
    int rsnx = root.get_nx();
    int rny  = root.get_ny();

    field2 R(rsnx, rny);

    for (int i = 0; i < csn; i++)
        buf_csn[i] = root(i, rny - 1);
    MPI_Allgatherv(buf_csn, csn, MPI_DOUBLE, buf_cn, cs_lengths, cs_displacements, MPI_DOUBLE, communicator);

    OPENMP_PARALLEL_FOR()
    for (int i = 0; i < csn; i++)
    {
        R(i, rny - 1) = 0.;
        for (int j = 0; j < cn; j++)
            R(i, rny - 1) += buf_cn[j] * value(i, j);
    }
    return R;
}

void SchurMat2DSlabX_down::construct(DomainSolver2D* branch_solver)
{
    field2 t_a(bsnx, bny);
    for (int i = 0; i < cn; i++)
    {
        t_a.clear();
        if (i >= cs_displacements[mpi_rank])
        {
            if (mpi_rank < mpi_size - 1 && i >= cs_displacements[mpi_rank + 1])
                continue;
            t_a(i - cs_displacements[mpi_rank], bny - 1) = 1.;
        }
        branch_solver->solve(t_a);
        for (int j = 0; j < csn; j++)
            value(j, i) = t_a(j, bny - 1);
    }
}

field2 SchurMat2DSlabX_down::operator*(const field2& root)
{
    int rsnx = root.get_nx();
    int rny  = root.get_ny();

    field2 R(rsnx, rny);

    for (int i = 0; i < csn; i++)
        buf_csn[i] = root(i, 0);
    MPI_Allgatherv(buf_csn, csn, MPI_DOUBLE, buf_cn, cs_lengths, cs_displacements, MPI_DOUBLE, communicator);

    OPENMP_PARALLEL_FOR()
    for (int i = 0; i < csn; i++)
    {
        R(i, 0) = 0.;
        for (int j = 0; j < cn; j++)
            R(i, 0) += buf_cn[j] * value(i, j);
    }
    return R;
}
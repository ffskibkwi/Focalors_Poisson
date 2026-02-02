#pragma once

#include "base/parallel/mpi/mpi_misc.h"
#include "base/pch.h"

#include <mpi.h>

class DomainSolver2D;

class SchurMat2DSlabX
{
protected:
    int    bnx; // branch domain nx
    int    bny; // branch domain ny
    int    cn;  // concat interface total n
    field2 value;

    // slab
    int      bsnx; // branch domain with slab decomposition nx
    int      csn;  // concat interface with slab decomposition total n
    MPI_Comm communicator;
    int      mpi_rank, mpi_size;
    int*     cs_lengths       = nullptr; // concat interface with slab decomposition's length for each rank
    int*     cs_displacements = nullptr; // concat interface with slab decomposition's displacement for each rank
    double*  buf_csn          = nullptr; // buffer size = csn
    double*  buf_cn           = nullptr; // buffer size = cn

    LocationType loc;

public:
    SchurMat2DSlabX(const Domain2DUniform* domain, const int _cn, MPI_Comm _communicator)
        : bnx(domain->get_nx())
        , bny(domain->get_ny())
        , cn(_cn)
        , communicator(_communicator)
    {
        MPI_Comm_rank(communicator, &mpi_rank);
        MPI_Comm_size(communicator, &mpi_size);

        bsnx = MPIUtils::get_slab_length(bnx, mpi_rank, mpi_size);
        csn  = MPIUtils::get_slab_length(cn, mpi_rank, mpi_size);

        value.init(csn, cn);

        buf_csn = new double[csn];
        buf_cn  = new double[cn];

        cs_lengths       = new int[mpi_size];
        cs_displacements = new int[mpi_size];
        for (int i = 0; i < mpi_size; i++)
        {
            cs_lengths[i]       = MPIUtils::get_slab_length(cn, i, mpi_size);
            cs_displacements[i] = MPIUtils::get_slab_displacement(cn, i, mpi_size);
        }
    }

    ~SchurMat2DSlabX()
    {
        delete[] buf_csn;
        delete[] buf_cn;

        delete[] cs_lengths;
        delete[] cs_displacements;
    }

    int          get_size() const { return cn; }
    LocationType get_loc() const { return loc; }
    void         print() { value.print(); }

    // Input branch solver should solve at slab x
    virtual void construct(DomainSolver2D* branch_solver) = 0;
    // Input root should alloc at slab x
    virtual field2 operator*(const field2& root) = 0;

    void set_name(const std::string& name) { value.set_name(name); }
    void write_csv(const std::string& directory);

    // Migrate from src to dest
    // The src pointer is only valid in src communicator
    // dest is only valid in dest communicator
    // Caller should ensure all processes in both comm will call this function
    friend void migrate_from(SchurMat2DSlabX* src, SchurMat2DSlabX* dest);
};

class SchurMat2DSlabX_left : public SchurMat2DSlabX
{
public:
    SchurMat2DSlabX_left(const Domain2DUniform* domain, MPI_Comm _communicator)
        : SchurMat2DSlabX(domain, domain->get_ny(), _communicator)
    {
        loc = LocationType::Left;
    }
    void   construct(DomainSolver2D* branch_solver) override;
    field2 operator*(const field2& root) override;
};

class SchurMat2DSlabX_right : public SchurMat2DSlabX
{
public:
    SchurMat2DSlabX_right(const Domain2DUniform* domain, MPI_Comm _communicator)
        : SchurMat2DSlabX(domain, domain->get_ny(), _communicator)
    {
        loc = LocationType::Right;
    }
    void   construct(DomainSolver2D* branch_solver) override;
    field2 operator*(const field2& root) override;
};

class SchurMat2DSlabX_up : public SchurMat2DSlabX
{
public:
    SchurMat2DSlabX_up(const Domain2DUniform* domain, MPI_Comm _communicator)
        : SchurMat2DSlabX(domain, domain->get_nx(), _communicator)
    {
        loc = LocationType::Up;
    }
    void   construct(DomainSolver2D* branch_solver) override;
    field2 operator*(const field2& root) override;
};

class SchurMat2DSlabX_down : public SchurMat2DSlabX
{
public:
    SchurMat2DSlabX_down(const Domain2DUniform* domain, MPI_Comm _communicator)
        : SchurMat2DSlabX(domain, domain->get_nx(), _communicator)
    {
        loc = LocationType::Down;
    }
    void   construct(DomainSolver2D* branch_solver) override;
    field2 operator*(const field2& root) override;
};

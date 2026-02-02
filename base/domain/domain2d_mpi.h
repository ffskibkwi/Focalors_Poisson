#pragma once

#include "domain2d.h"

#include <mpi.h>

class Domain2DMPIUniform : public Domain2DUniform
{
public:
    Domain2DMPIUniform(MPI_Comm _communicator = MPI_COMM_WORLD);
    Domain2DMPIUniform(const std::string& in_name, MPI_Comm _communicator = MPI_COMM_WORLD);
    Domain2DMPIUniform(int                in_nx,
                       int                in_ny,
                       double             in_lx,
                       double             in_ly,
                       const std::string& in_name,
                       MPI_Comm           _communicator = MPI_COMM_WORLD);
    Domain2DMPIUniform(int in_nx, int in_ny, const std::string& in_name, MPI_Comm _communicator = MPI_COMM_WORLD);

    int get_uuid() { return uuid; }

private:
    void sync_uuid();

    // Domains point to same physical position in different process,
    // are different address in process local memory.
    // So we need a uuid to track domain.
    int uuid = -1;

    // slab
    MPI_Comm communicator = MPI_COMM_WORLD;
    int      mpi_rank, mpi_size;

    static int s_counter;
};
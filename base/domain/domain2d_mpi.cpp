#include "domain2d_mpi.h"

int Domain2DMPIUniform::s_counter = 0;

Domain2DMPIUniform::Domain2DMPIUniform(MPI_Comm _communicator)
    : communicator(_communicator)
{
    MPI_Comm_rank(communicator, &mpi_rank);
    MPI_Comm_size(communicator, &mpi_size);

    uuid = s_counter++;
}

Domain2DMPIUniform::Domain2DMPIUniform(const std::string& in_name, MPI_Comm _communicator)
    : Domain2DUniform(in_name)
    , communicator(_communicator)
{
    sync_uuid();
}

Domain2DMPIUniform::Domain2DMPIUniform(int                in_nx,
                                       int                in_ny,
                                       double             in_lx,
                                       double             in_ly,
                                       const std::string& in_name,
                                       MPI_Comm           _communicator)
    : Domain2DUniform(in_nx, in_ny, in_lx, in_ly, in_name)
    , communicator(_communicator)
{
    sync_uuid();
}

Domain2DMPIUniform::Domain2DMPIUniform(int in_nx, int in_ny, const std::string& in_name, MPI_Comm _communicator)
    : Domain2DUniform(in_nx, in_ny, in_name)
    , communicator(_communicator)
{
    sync_uuid();
}

void Domain2DMPIUniform::sync_uuid()
{
    MPI_Comm_rank(communicator, &mpi_rank);
    MPI_Comm_size(communicator, &mpi_size);

    if (mpi_rank == 0)
        uuid = s_counter++;

    MPI_Bcast(&uuid, 1, MPI_INT, 0, communicator);
}
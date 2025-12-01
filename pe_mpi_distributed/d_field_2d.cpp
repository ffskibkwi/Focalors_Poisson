#include "d_field_2d.h"
#include <iostream>

DField2D::DField2D(DDomain2D* in_domain, std::string in_name)
    : domain(in_domain)
    , name(in_name)
{
    if (domain)
    {
        local_data.init(domain->get_local_nx(), domain->get_local_ny(), name + "_local");
    }
}

void DField2D::zeros() { local_data.clear(); }

void DField2D::get_global_coord(int i_local, int j_local, int& i_global, int& j_global) const
{
    if (!domain)
        return;
    i_global = domain->get_local_i_start() + i_local;
    j_global = j_local; // For i-slab, j is not split
}

void DField2D::scatter_from_global(const field2& global_field)
{
    if (!domain)
        return;

    int      rank = domain->get_rank();
    int      size = domain->get_size();
    MPI_Comm comm = domain->get_comm();
    int      ny   = domain->get_global_ny();

    const std::vector<int>& i_counts = domain->get_i_counts();
    const std::vector<int>& i_displs = domain->get_i_displs();

    std::vector<int> sendcounts(size), displs(size);
    for (int r = 0; r < size; ++r)
    {
        sendcounts[r] = i_counts[r] * ny;
        displs[r]     = i_displs[r] * ny;
    }

    // Note: MPI_Scatterv expects sendbuf to be valid on root
    const double* sendbuf = (rank == 0) ? global_field.value : nullptr;

    MPI_Scatterv(
        sendbuf, sendcounts.data(), displs.data(), MPI_DOUBLE, local_data.value, sendcounts[rank], MPI_DOUBLE, 0, comm);
}

void DField2D::gather_to_global(field2& global_field) const
{
    if (!domain)
        return;

    int      rank = domain->get_rank();
    int      size = domain->get_size();
    MPI_Comm comm = domain->get_comm();
    int      nx   = domain->get_global_nx();
    int      ny   = domain->get_global_ny();

    const std::vector<int>& i_counts = domain->get_i_counts();
    const std::vector<int>& i_displs = domain->get_i_displs();

    std::vector<int> recvcounts(size), displs(size);
    for (int r = 0; r < size; ++r)
    {
        recvcounts[r] = i_counts[r] * ny;
        displs[r]     = i_displs[r] * ny;
    }

    if (rank == 0)
    {
        global_field.init(nx, ny);
    }

    double* recvbuf = (rank == 0) ? global_field.value : nullptr;

    MPI_Gatherv(
        local_data.value, recvcounts[rank], MPI_DOUBLE, recvbuf, recvcounts.data(), displs.data(), MPI_DOUBLE, 0, comm);
}

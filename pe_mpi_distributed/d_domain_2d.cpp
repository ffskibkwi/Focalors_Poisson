#include "d_domain_2d.h"

DDomain2D::DDomain2D(MPI_Comm in_comm, int in_global_nx, int in_global_ny, double in_hx, double in_hy)
    : comm(in_comm), global_nx(in_global_nx), global_ny(in_global_ny), hx(in_hx), hy(in_hy)
{
    if (comm != MPI_COMM_NULL)
    {
        MPI_Comm_rank(comm, &rank);
        MPI_Comm_size(comm, &size);
    }
    else
    {
        rank = -1;
        size = 0;
    }

    if (size > 0)
    {
        split_1d(global_nx, size, i_counts, i_displs);
        local_nx = i_counts[rank];
        local_ny = global_ny; // Full Y range for i-slab
    }
}

void DDomain2D::get_j_decomposition(std::vector<int>& out_j_counts, std::vector<int>& out_j_displs) const
{
    if (size > 0)
    {
        split_1d(global_ny, size, out_j_counts, out_j_displs);
    }
}

void DDomain2D::split_1d(int n, int p, std::vector<int>& counts, std::vector<int>& displs)
{
    counts.resize(p);
    displs.resize(p);
    int base = n / p;
    int rem  = n % p;
    int off  = 0;
    for (int r = 0; r < p; ++r)
    {
        counts[r] = base + (r < rem ? 1 : 0);
        displs[r] = off;
        off += counts[r];
    }
}

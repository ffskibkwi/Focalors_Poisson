#include "distribute_slab.h"
#include "mpi_misc.h"
#include <algorithm>
#include <vector>

namespace MPIUtils
{
    void redistribute_2d_slab_sync(field2* matrix_src, field2* matrix_dest, MPI_Comm comm_src, MPI_Comm comm_dest)
    {
        int world_rank, world_size;
        MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
        MPI_Comm_size(MPI_COMM_WORLD, &world_size);

        std::vector<int> src_ranks_in_world(world_size, -1);
        std::vector<int> dest_ranks_in_world(world_size, -1);

        int my_rank_src = -1, size_src = 0;
        int my_rank_dest = -1, size_dest = 0;
        int local_nx = 0, ny_full = 0;

        // Check source pointer and communicator
        if (comm_src != MPI_COMM_NULL && matrix_src != nullptr)
        {
            MPI_Comm_rank(comm_src, &my_rank_src);
            MPI_Comm_size(comm_src, &size_src);
            local_nx = matrix_src->get_nx();
            ny_full  = matrix_src->get_ny();
        }

        // Check destination pointer and communicator
        if (comm_dest != MPI_COMM_NULL && matrix_dest != nullptr)
        {
            MPI_Comm_rank(comm_dest, &my_rank_dest);
            MPI_Comm_size(comm_dest, &size_dest);
            ny_full = matrix_dest->get_ny();
        }

        // Synchronize metadata across MPI_COMM_WORLD
        MPI_Allgather(&my_rank_src, 1, MPI_INT, src_ranks_in_world.data(), 1, MPI_INT, MPI_COMM_WORLD);
        MPI_Allgather(&my_rank_dest, 1, MPI_INT, dest_ranks_in_world.data(), 1, MPI_INT, MPI_COMM_WORLD);

        int nx_full = 0;
        MPI_Allreduce(&local_nx, &nx_full, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(MPI_IN_PLACE, &ny_full, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
        MPI_Allreduce(MPI_IN_PLACE, &size_src, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
        MPI_Allreduce(MPI_IN_PLACE, &size_dest, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);

        std::vector<int> src_to_world(size_src);
        std::vector<int> dest_to_world(size_dest);
        for (int i = 0; i < world_size; ++i)
        {
            if (src_ranks_in_world[i] != -1)
                src_to_world[src_ranks_in_world[i]] = i;
            if (dest_ranks_in_world[i] != -1)
                dest_to_world[dest_ranks_in_world[i]] = i;
        }

        std::vector<MPI_Request> requests;

        // Sender Side logic
        if (comm_src != MPI_COMM_NULL && matrix_src != nullptr)
        {
            int my_src_start = get_slab_displacement(nx_full, my_rank_src, size_src);
            int my_src_nx    = matrix_src->get_nx();
            int my_src_end   = my_src_start + my_src_nx;

            for (int p_dest = 0; p_dest < size_dest; ++p_dest)
            {
                int target_start = get_slab_displacement(nx_full, p_dest, size_dest);
                int target_nx    = get_slab_length(nx_full, p_dest, size_dest);
                int target_end   = target_start + target_nx;

                int overlap_start = std::max(my_src_start, target_start);
                int overlap_end   = std::min(my_src_end, target_end);

                if (overlap_start < overlap_end)
                {
                    int         world_target = dest_to_world[p_dest];
                    double*     send_ptr     = matrix_src->get_ptr(overlap_start - my_src_start, 0);
                    MPI_Request req;
                    MPI_Isend(send_ptr,
                              (overlap_end - overlap_start) * ny_full,
                              MPI_DOUBLE,
                              world_target,
                              101,
                              MPI_COMM_WORLD,
                              &req);
                    requests.push_back(req);
                }
            }
        }

        // Receiver Side logic
        if (comm_dest != MPI_COMM_NULL && matrix_dest != nullptr)
        {
            int my_dest_start = get_slab_displacement(nx_full, my_rank_dest, size_dest);
            int my_dest_nx    = matrix_dest->get_nx();
            int my_dest_end   = my_dest_start + my_dest_nx;

            for (int p_src = 0; p_src < size_src; ++p_src)
            {
                int src_start = get_slab_displacement(nx_full, p_src, size_src);
                int src_nx    = get_slab_length(nx_full, p_src, size_src);
                int src_end   = src_start + src_nx;

                int overlap_start = std::max(my_dest_start, src_start);
                int overlap_end   = std::min(my_dest_end, src_end);

                if (overlap_start < overlap_end)
                {
                    int         world_src = src_to_world[p_src];
                    double*     recv_ptr  = matrix_dest->get_ptr(overlap_start - my_dest_start, 0);
                    MPI_Request req;
                    MPI_Irecv(recv_ptr,
                              (overlap_end - overlap_start) * ny_full,
                              MPI_DOUBLE,
                              world_src,
                              101,
                              MPI_COMM_WORLD,
                              &req);
                    requests.push_back(req);
                }
            }
        }

        if (!requests.empty())
            MPI_Waitall(static_cast<int>(requests.size()), requests.data(), MPI_STATUSES_IGNORE);
    }
} // namespace MPIUtils
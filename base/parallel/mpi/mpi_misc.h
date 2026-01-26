#pragma once

#include <cstdlib> // for std::getenv
#include <iostream>
#include <string>

namespace MPIUtils
{
    bool is_mpi();

    template<typename Func>
    void do_once(bool is_mpi, int mpi_rank, Func&& func)
    {
        if (is_mpi)
        {
            if (mpi_rank == 0)
            {
                func();
            }
        }
        else
        {
            func();
        }
    }

    void bcast_string(std::string& src_recv_str, int char_buf_size = 100);

    int get_slab_length(int total_length, int mpi_rank, int mpi_size);
    int get_slab_displacement(int total_length, int mpi_rank, int mpi_size);
} // namespace MPIUtils
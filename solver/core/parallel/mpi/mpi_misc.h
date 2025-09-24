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
} // namespace MPIUtils
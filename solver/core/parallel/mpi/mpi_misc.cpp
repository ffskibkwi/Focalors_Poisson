#include "mpi_misc.h"

#include <mpi.h>

namespace MPIUtils
{
    bool is_mpi()
    {
        static const char* env_vars[] = {
            "OMPI_COMM_WORLD_SIZE", // Open MPI
            "PMI_RANK",             // SLURM and others
            "MPICH_RANK",           // MPICH
            "PMI_SIZE"              // General MPI
        };

        for (const char* var : env_vars)
        {
            if (std::getenv(var) != nullptr)
            {
                return true;
            }
        }

        return false;
    }

    void bcast_string(std::string& src_recv_str, int char_buf_size)
    {
        int mpi_rank = 0;

        MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);

        char* char_buf = new char[char_buf_size];
        for (int i = 0; i < char_buf_size; ++i)
        {
            char_buf[i] = '\0';
        }
        if (mpi_rank == 0)
        {
            for (int i = 0; i < src_recv_str.size(); ++i)
            {
                char_buf[i] = src_recv_str[i];
            }
        }

        MPI_Bcast(char_buf, char_buf_size, MPI_CHAR, 0, MPI_COMM_WORLD);
        src_recv_str = std::string(char_buf);

        delete[] char_buf;
    }
} // namespace MPIUtils
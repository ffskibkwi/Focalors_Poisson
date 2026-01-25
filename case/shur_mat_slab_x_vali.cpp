#include "pe/concat/domain_solver.h"
#include "pe/concat/schur_mat2d_slab_x.h"

#include <unistd.h>

class DomainSolver2DTest : public DomainSolver2D
{
public:
    int disp = 0;

    DomainSolver2DTest(int _disp)
        : disp(_disp)
    {}
    ~DomainSolver2DTest() {}
    void solve(field2& f, bool is_debugmode = true)
    {
        for (int i = 0; i < f.get_nx(); i++)
        {
            for (int j = 0; j < f.get_ny(); j++)
            {
                f(i, j) *= (i + disp) * f.get_ny() + j + 1;
            }
        }
    }
    bool is_comm_root() const { return true; }

    double get_hx() const { return 1; };
    double get_hy() const { return 1; };
};

int main(int argc, char* argv[])
{
    MPI_Init(&argc, &argv);

    // debug
    volatile int ii = 0;
    while (ii == 0)
        sleep(1);

    int    nx          = 4;
    int    ny          = 3;
    int    neighbor_nx = 6;
    int    neighbor_ny = 3;
    double lx          = 1.0;
    double ly          = 1.0;

    int mpi_rank = 0, mpi_size = 1;

    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);

    int nx_slab = (mpi_rank == mpi_size - 1) ? (nx - nx / mpi_size * mpi_rank) : nx / mpi_size;
    int nx_disp = nx / mpi_size * mpi_rank;

    Domain2DUniform neighbor_domain(neighbor_nx, neighbor_ny, lx, ly, "Test");

    field2 f(nx_slab, ny);

    for (int i = 0; i < f.get_nx(); i++)
    {
        for (int j = 0; j < f.get_ny(); j++)
        {
            f(i, j) = (i + nx_disp) * f.get_ny() + j + 1;
        }
    }

    // for (int i = 0; i < mpi_size; i++)
    // {
    //     if (i == mpi_rank)
    //     {
    //         std::cout << "rank " << mpi_rank << std::endl;
    //         f.print();
    //     }
    //     MPI_Barrier(MPI_COMM_WORLD);
    // }

    // if (mpi_rank == 0)
    //     std::cout << "---------------" << std::endl;

    SchurMat2DSlabX_left shur(&neighbor_domain, MPI_COMM_WORLD);
    DomainSolver2DTest   solver(nx_disp);
    shur.construct(&solver);
    f = shur * f;

    // for (int i = 0; i < mpi_size; i++)
    // {
    //     if (i == mpi_rank)
    //     {
    //         std::cout << "rank " << mpi_rank << std::endl;
    //         shur.print();
    //     }
    //     MPI_Barrier(MPI_COMM_WORLD);
    // }

    // if (mpi_rank == 0)
    //     std::cout << "---------------" << std::endl;

    // for (int i = 0; i < mpi_size; i++)
    // {
    //     if (i == mpi_rank)
    //     {
    //         std::cout << "rank " << mpi_rank << std::endl;
    //         f.print();
    //     }
    //     MPI_Barrier(MPI_COMM_WORLD);
    // }

    MPI_Finalize();

    return 0;
}
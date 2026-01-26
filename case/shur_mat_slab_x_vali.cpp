#include "base/parallel/mpi/mpi_misc.h"
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
    void solve(field2& f)
    {
        for (int i = 0; i < f.get_nx(); i++)
        {
            for (int j = 0; j < f.get_ny(); j++)
            {
                f(i, j) *= (i + disp) * f.get_ny() + j + 1;
            }
        }
    }
};

template<typename SchurType>
void test_schur_direction(const std::string& label, Domain2DUniform& neighbor_domain, int nx, int ny)
{
    int mpi_rank = 0, mpi_size = 1;

    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);

    if (mpi_rank == 0)
        std::cout << "\n========== Testing " << label << " ==========" << std::endl;

    int neighbor_nx      = neighbor_domain.get_nx();
    int nx_slab          = MPIUtils::get_slab_length(nx, mpi_rank, mpi_size);
    int nx_disp          = MPIUtils::get_slab_displacement(nx, mpi_rank, mpi_size);
    int neighbor_nx_disp = MPIUtils::get_slab_displacement(neighbor_nx, mpi_rank, mpi_size);

    field2 f(nx_slab, ny);

    for (int i = 0; i < f.get_nx(); i++)
    {
        for (int j = 0; j < f.get_ny(); j++)
        {
            f(i, j) = (i + nx_disp) * f.get_ny() + j + 1;
        }
    }

    if (mpi_rank == 0)
        std::cout << "Original Field:" << std::endl;
    for (int i = 0; i < mpi_size; i++)
    {
        if (i == mpi_rank)
        {
            std::cout << "rank " << mpi_rank << std::endl;
            f.print();
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }

    SchurType          shur(&neighbor_domain, MPI_COMM_WORLD);
    DomainSolver2DTest solver(neighbor_nx_disp);
    shur.construct(&solver);

    f = shur * f;

    if (mpi_rank == 0)
        std::cout << "Schur Matrix Structure:" << std::endl;
    for (int i = 0; i < mpi_size; i++)
    {
        if (i == mpi_rank)
        {
            std::cout << "rank " << mpi_rank << std::endl;
            shur.print();
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }

    if (mpi_rank == 0)
        std::cout << "Result Field (after Schur * f):" << std::endl;
    for (int i = 0; i < mpi_size; i++)
    {
        if (i == mpi_rank)
        {
            std::cout << "rank " << mpi_rank << std::endl;
            f.print();
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }
}

int main(int argc, char* argv[])
{
    MPI_Init(&argc, &argv);

    int    nx          = 4;
    int    ny          = 3;
    int    neighbor_nx = 6;
    int    neighbor_ny = 3;
    double lx          = 1.0;
    double ly          = 1.0;

    Domain2DUniform neighbor_domain_x(neighbor_nx, neighbor_ny, lx, ly, "TestDomain");
    Domain2DUniform neighbor_domain_y(neighbor_ny, neighbor_nx, ly, lx, "TestDomain");

    test_schur_direction<SchurMat2DSlabX_left>("LEFT", neighbor_domain_x, nx, ny);
    test_schur_direction<SchurMat2DSlabX_right>("RIGHT", neighbor_domain_x, nx, ny);
    test_schur_direction<SchurMat2DSlabX_up>("UP", neighbor_domain_y, ny, nx);
    test_schur_direction<SchurMat2DSlabX_down>("DOWN", neighbor_domain_y, ny, nx);

    MPI_Finalize();

    return 0;
}
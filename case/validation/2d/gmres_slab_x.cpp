#include "base/domain/domain2d.h"
#include "base/domain/geometry2d.h"
#include "base/domain/geometry_tree.hpp"
#include "base/domain/variable2d.h"
#include "base/field/field2.h"
#include "base/parallel/mpi/mpi_misc.h"
#include "pe/concat/gmres_solver2d_slab_x.h"

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

void test(LocationType neighbor_loc)
{
    int mpi_rank = 0, mpi_size = 1;

    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);

    if (mpi_rank == 0)
        std::cout << "\n========== Testing " << neighbor_loc << " ==========" << std::endl;

    int nx          = 3;
    int ny          = 5;
    int neighbor_n  = 4;
    int neighbor_nx = nx, neighbor_ny = ny;
    if (neighbor_loc == LocationType::Left || neighbor_loc == LocationType::Right)
        neighbor_nx = neighbor_n;
    else if (neighbor_loc == LocationType::Down || neighbor_loc == LocationType::Up)
        neighbor_ny = neighbor_n;

    int    m       = 2;
    double tol     = 1e-3;
    int    maxIter = 2;

    Geometry2D geo_tee;

    int nx_slab          = MPIUtils::get_slab_length(nx, mpi_rank, mpi_size);
    int nx_disp          = MPIUtils::get_slab_displacement(nx, mpi_rank, mpi_size);
    int neighbor_nx_disp = MPIUtils::get_slab_displacement(neighbor_nx, mpi_rank, mpi_size);

    Domain2DUniform root(nx, ny, "root");
    Domain2DUniform neighbor(neighbor_nx, neighbor_ny, "neighbor");
    field2          p_root(nx_slab, ny, "p_root");

    GMRESSolver2DSlabX  gmres(&root, m, tol, maxIter);
    DomainSolver2DTest* test_solver_root     = new DomainSolver2DTest(nx_disp);
    DomainSolver2DTest* test_solver_neighbor = new DomainSolver2DTest(neighbor_nx_disp);
    gmres.set_solver(test_solver_root);

    std::unordered_map<LocationType, Domain2DUniform*> adj;
    adj[neighbor_loc] = &neighbor;

    std::unordered_map<Domain2DUniform*, DomainSolver2D*> solver_map;
    solver_map[&root]     = &gmres;
    solver_map[&neighbor] = test_solver_neighbor;

    gmres.schur_mat_construct(adj, solver_map);

    for (int i = 0; i < p_root.get_nx(); i++)
    {
        for (int j = 0; j < p_root.get_ny(); j++)
        {
            p_root(i, j) = (i + nx_disp) * p_root.get_ny() + j + 1;
        }
    }

    gmres.solve(p_root);

    for (int i = 0; i < mpi_size; i++)
    {
        if (i == mpi_rank)
        {
            std::cout << "rank " << mpi_rank << std::endl;
            p_root.print();
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }

    delete test_solver_neighbor;
}

int main(int argc, char* argv[])
{
    MPI_Init(&argc, &argv);

    test(LocationType::Left);
    test(LocationType::Right);
    test(LocationType::Down);
    test(LocationType::Up);

    MPI_Finalize();

    return 0;
}
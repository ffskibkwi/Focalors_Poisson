#include "base/domain/variable2d_slab_x.h"

int main(int argc, char* argv[])
{
    MPI_Init(&argc, &argv);

    int mpi_rank = 0, mpi_size = 1;

    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);

    int    rank = 4;
    int    n1 = rank, n2 = 2 * rank, n3 = 4 * rank;
    int    m4 = 2 * rank, m2 = 4 * rank, m5 = rank;
    double H = 1.0 / rank;

    Domain2DUniform T1(n1, m2, "T1");
    Domain2DUniform T2(n2, m2, "T2");
    Domain2DUniform T3(n3, m2, "T3");
    Domain2DUniform T4(n2, m4, "T4");
    Domain2DUniform T5(n2, m5, "T5");

    Geometry2D geo;
    geo.connect(&T2, LocationType::Left, &T1);
    geo.connect(&T2, LocationType::Right, &T3);
    geo.connect(&T2, LocationType::Down, &T4);
    geo.connect(&T2, LocationType::Up, &T5);

    geo.set_global_spatial_step(H, H);

    geo.check();
    geo.solve_prepare();

    Variable2DSlabX p;
    p.set_geometry(geo);

    for (int level = 0; level < 2; level++)
    {
        if (mpi_rank == 0)
            std::cout << "----------------------------" << std::endl;

        for (int i = 0; i < mpi_size; i++)
        {
            if (i == mpi_rank)
            {
                std::cout << "rank " << mpi_rank << ' ';
                std::cout << p.hierarchical_slab_parents[level]->name << std::endl;
            }
            MPI_Barrier(MPI_COMM_WORLD);
        }
    }

    MPI_Finalize();

    return 0;
}
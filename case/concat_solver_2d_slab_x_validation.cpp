#include "base/domain/domain2d_mpi.h"
#include "base/domain/geometry2d.h"
#include "base/domain/geometry_tree.hpp"
#include "base/domain/variable2d.h"
#include "base/field/field2.h"
#include "base/location_boundary.h"
#include "base/parallel/mpi/mpi_misc.h"
#include "io/config.h"
#include "io/csv_writer_2d.h"
#include "pe/concat/concat_solver2d_slab_x.h"

#include <unistd.h> // for sleep

void fill(field2& f, int disp)
{
    for (int i = 0; i < f.get_nx(); i++)
    {
        for (int j = 0; j < f.get_ny(); j++)
        {
            f(i, j) = (i + disp) * f.get_ny() + j + 1;
        }
    }
}

void print_mpi(field2& f, int mpi_rank, int mpi_size)
{
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

    int mpi_rank = 0, mpi_size = 1;

    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);

    // debug
    // volatile int ii = 0; // use volatile to avoid compiler optimatize out variable
    // while (ii == 0)
    //     sleep(1);

    Geometry2D         geo_tee;
    EnvironmentConfig* env_config = new EnvironmentConfig();

    int nx_2 = 3;
    int ny_2 = 3;
    int nx_1 = 6;
    int ny_3 = 9;

    int nx_2_disp = MPIUtils::get_slab_displacement(nx_2, mpi_rank, mpi_size);
    int nx_1_disp = MPIUtils::get_slab_displacement(nx_2, mpi_rank, mpi_size);

    Domain2DUniformMPI T2(nx_2, ny_2, 1.0, 1.0, "T2");
    Domain2DUniformMPI T1("T1");
    T1.set_nx(nx_1);
    T1.set_lx(2.0);
    Domain2DUniformMPI T3("T3");
    T3.set_ny(ny_3);
    T3.set_ly(3.0);

    geo_tee.add_domain({&T1, &T2, &T3});

    geo_tee.connect(&T2, LocationType::Left, &T1);
    geo_tee.connect(&T2, LocationType::Down, &T3);

    geo_tee.check();
    geo_tee.solve_prepare();

    Variable2DSlabX v("v");
    v.set_geometry(geo_tee);
    field2 v_T1, v_T2, v_T3;
    v.set_center_field(&T1, v_T1);
    v.set_center_field(&T2, v_T2);
    v.set_center_field(&T3, v_T3);

    v.fill_boundary_type(PDEBoundaryType::Dirichlet);

    fill(v_T1, nx_1_disp);
    fill(v_T2, nx_2_disp);
    fill(v_T3, nx_2_disp);

    v.print_slab_info();

    ConcatPoissonSolver2DSlabX solver(&v, env_config);
    solver.solve();

    print_mpi(v_T1, mpi_rank, mpi_size);
    std::cout << "--------------------------" << std::endl;
    print_mpi(v_T2, mpi_rank, mpi_size);
    std::cout << "--------------------------" << std::endl;
    print_mpi(v_T3, mpi_rank, mpi_size);

    MPI_Finalize();

    return 0;
}
#include "base/config.h"
#include "base/domain/domain2d_mpi.h"
#include "base/domain/geometry2d.h"
#include "base/domain/geometry_tree.hpp"
#include "base/domain/variable2d.h"
#include "base/field/field2.h"
#include "base/location_boundary.h"
#include "base/parallel/mpi/mpi_misc.h"
#include "pe/concat/concat_solver2d_slab_x.h"

#include <unistd.h> // for sleep

void fill(field2& f, Variable2DSlabX* var, Domain2DUniformMPI* domain)
{
    int disp = 0;
    if (var->slab_parent_to_level.find(domain->get_uuid()) != var->slab_parent_to_level.end())
        disp = var->hierarchical_slab_disps[var->slab_parent_to_level[domain->get_uuid()]];

    for (int i = 0; i < f.get_nx(); i++)
    {
        for (int j = 0; j < f.get_ny(); j++)
        {
            f(i, j) = (i + disp) * f.get_ny() + j + 1;
        }
    }
}

/**
 *
 * y
 * ▲
 * │
 * │      ┌──────┐
 * │      │      │
 * │      │  A5  │
 * │      │      │
 * ├──────┼──────┼──────┬──────┐
 * │      │      │      │      │
 * │  A2  │  A1  │  A3  │  A6  │
 * │      │      │      │      │
 * ├──────┼──────┼──────┴──────┘
 * │      │      │
 * │      │  A4  │
 * │      │      │
 * └──────┴──────┴──────────►
 * O                         x
 *
 */
int main(int argc, char* argv[])
{
    MPI_Init(&argc, &argv);

    int mpi_rank = 0, mpi_size = 1;

    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);

    // // debug
    // volatile int ii = 0; // use volatile to avoid compiler optimatize out variable
    // while (ii == 0)
    //     sleep(1);

    Geometry2D geo;

    int nx_1 = 8;
    int ny_1 = 7;
    int nx_2 = 7;
    int nx_3 = 8;
    int ny_4 = 8;
    int ny_5 = 7;
    int nx_6 = 8;

    Domain2DUniformMPI T1(nx_1, ny_1, nx_1, ny_1, "T1");
    Domain2DUniformMPI T2("T2");
    T2.set_nx(nx_2);
    T2.set_lx(nx_2);
    Domain2DUniformMPI T3("T3");
    T3.set_nx(nx_3);
    T3.set_lx(nx_3);
    Domain2DUniformMPI T4("T4");
    T4.set_ny(ny_4);
    T4.set_ly(ny_4);
    Domain2DUniformMPI T5("T5");
    T5.set_ny(ny_5);
    T5.set_ly(ny_5);
    Domain2DUniformMPI T6("T6");
    T6.set_nx(nx_6);
    T6.set_lx(nx_6);

    geo.add_domain({&T1, &T2, &T3, &T4, &T5, &T6});

    geo.connect(&T1, LocationType::Left, &T2);
    geo.connect(&T1, LocationType::Right, &T3);
    geo.connect(&T1, LocationType::Down, &T4);
    geo.connect(&T1, LocationType::Up, &T5);
    geo.connect(&T3, LocationType::Right, &T6);

    geo.check();
    geo.solve_prepare();

    Variable2DSlabX v("v");
    v.set_geometry(geo);
    field2 v_T1, v_T2, v_T3, v_T4, v_T5, v_T6;
    v.set_center_field(&T1, v_T1);
    v.set_center_field(&T2, v_T2);
    v.set_center_field(&T3, v_T3);
    v.set_center_field(&T4, v_T4);
    v.set_center_field(&T5, v_T5);
    v.set_center_field(&T6, v_T6);

    v.fill_boundary_type(PDEBoundaryType::Dirichlet);

    for (auto kv : v.field_map)
        fill(*kv.second, &v, static_cast<Domain2DUniformMPI*>(kv.first));

    auto print_field = [&](field2& f) {
        std::cout << "--------------------------" << std::endl;
        f.print();
    };

    std::array<field2*, 6> v_to_print = {&v_T1, &v_T2, &v_T3, &v_T4, &v_T5, &v_T6};
    for (auto f : v_to_print)
        print_field(*f);

    std::cout << "-----------solve---------------" << std::endl;

    ConcatPoissonSolver2DSlabX solver(&v);
    solver.solve();

    for (auto f : v_to_print)
        print_field(*f);

    MPI_Finalize();

    return 0;
}
#include "base/domain/domain2d.h"
#include "base/domain/geometry2d.h"
#include "base/domain/geometry_tree.hpp"
#include "base/domain/variable2d.h"
#include "base/field/field2.h"
#include "base/location_boundary.h"
#include "io/config.h"
#include "io/csv_writer_2d.h"
#include "pe/concat/concat_solver2d.h"

void fill(field2& f)
{
    for (int i = 0; i < f.get_nx(); i++)
    {
        for (int j = 0; j < f.get_ny(); j++)
        {
            f(i, j) = i * f.get_ny() + j + 1;
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
    Geometry2D         geo_tee;
    EnvironmentConfig* env_config = new EnvironmentConfig();

    int nx_1 = 8;
    int ny_1 = 7;
    int nx_2 = 7;
    int nx_3 = 8;
    int ny_4 = 8;
    int ny_5 = 7;
    int nx_6 = 8;

    Domain2DUniform T1(nx_1, ny_1, nx_1, ny_1, "T1");
    Domain2DUniform T2("T2");
    T2.set_nx(nx_2);
    T2.set_lx(nx_2);
    Domain2DUniform T3("T3");
    T3.set_nx(nx_3);
    T3.set_lx(nx_3);
    Domain2DUniform T4("T4");
    T4.set_ny(ny_4);
    T4.set_ly(ny_4);
    Domain2DUniform T5("T5");
    T5.set_ny(ny_5);
    T5.set_ly(ny_5);
    Domain2DUniform T6("T6");
    T6.set_nx(nx_6);
    T6.set_lx(nx_6);

    geo_tee.add_domain({&T1, &T2, &T3, &T4, &T5, &T6});

    geo_tee.connect(&T1, LocationType::Left, &T2);
    geo_tee.connect(&T1, LocationType::Right, &T3);
    geo_tee.connect(&T1, LocationType::Down, &T4);
    geo_tee.connect(&T1, LocationType::Up, &T5);
    geo_tee.connect(&T3, LocationType::Right, &T6);

    geo_tee.check();
    geo_tee.solve_prepare();

    Variable2D v("v");
    v.set_geometry(geo_tee);
    field2 v_T1, v_T2, v_T3, v_T4, v_T5, v_T6;
    v.set_center_field(&T1, v_T1);
    v.set_center_field(&T2, v_T2);
    v.set_center_field(&T3, v_T3);
    v.set_center_field(&T4, v_T4);
    v.set_center_field(&T5, v_T5);
    v.set_center_field(&T6, v_T6);

    v.fill_boundary_type(PDEBoundaryType::Dirichlet);

    for (auto kv : v.field_map)
        fill(*kv.second);

    auto print_field = [&](field2& f) {
        std::cout << "--------------------------" << std::endl;
        f.print();
    };

    std::array<field2*, 6> v_to_print = {&v_T1, &v_T2, &v_T3, &v_T4, &v_T5, &v_T6};
    for (auto f : v_to_print)
        print_field(*f);

    std::cout << "-----------solve---------------" << std::endl;

    ConcatPoissonSolver2D solver(&v, env_config);
    solver.solve();

    for (auto f : v_to_print)
        print_field(*f);

    return 0;
}
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

int main(int argc, char* argv[])
{
    Geometry2D         geo_tee;
    EnvironmentConfig* env_config = new EnvironmentConfig();
    env_config->debug_concat      = true;

    int nx_2 = 3;
    int ny_2 = 3;
    int nx_1 = 6;
    int ny_3 = 9;

    Domain2DUniform T2(nx_2, ny_2, 1.0, 1.0, "T2");
    Domain2DUniform T1("T1");
    T1.set_nx(nx_1);
    T1.set_lx(2.0);
    Domain2DUniform T3("T3");
    T3.set_ny(ny_3);
    T3.set_ly(3.0);

    geo_tee.add_domain({&T1, &T2, &T3});

    geo_tee.connect(&T2, LocationType::Left, &T1);
    geo_tee.connect(&T2, LocationType::Down, &T3);

    geo_tee.check();
    geo_tee.solve_prepare();

    Variable2D v("v");
    v.set_geometry(geo_tee);
    field2 v_T1, v_T2, v_T3;
    v.set_center_field(&T1, v_T1);
    v.set_center_field(&T2, v_T2);
    v.set_center_field(&T3, v_T3);

    v.fill_boundary_type(PDEBoundaryType::Dirichlet);

    fill(v_T1);
    fill(v_T2);
    fill(v_T3);

    auto print_field = [&](field2& f) {
        std::cout << "--------------------------" << std::endl;
        std::cout << f.get_name() << std::endl;
        f.print();
    };
    print_field(v_T1);
    print_field(v_T2);
    print_field(v_T3);

    ConcatPoissonSolver2D solver(&v, env_config);
    solver.solve();

    std::cout << "--------------------------" << std::endl;
    v_T1.print();
    std::cout << "--------------------------" << std::endl;
    v_T2.print();
    std::cout << "--------------------------" << std::endl;
    v_T3.print();

    return 0;
}
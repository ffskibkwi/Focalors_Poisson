#include "core/domain/geometry2d.h"
#include "core/domain/domain2d.h"
#include "core/domain/variable.h"
#include "core/base/field2.h"

#include "core/boundary/boundary_type.h"

#include "pe/concat_new/concat_solver2d.h"

int main(int argc, char* argv[])
{

    /*
            C
            |
        B---A
    */
    int nx_A = 10;
    int ny_A = 10;
    double lx_A = 1.0;
    double ly_A = 1.0;
    int nx_B = 10;
    int ny_B = 10;
    double lx_B = 1.0;
    double ly_B = 1.0;
    int nx_C = 10;
    int ny_C = 10;
    double lx_C = 1.0;
    double ly_C = 1.0;

    Geometry2D geo;
    Domain2DUniform A(nx_A, ny_A, lx_A, ly_A, "A");
    Domain2DUniform B("B");
    Domain2DUniform C("C");

    C.set_ly(ly_C);
    C.set_ny(ny_C);

    B.set_lx(lx_B);
    B.set_nx(nx_B);

    geo.add_domain(A);
    geo.add_domain(B);
    geo.add_domain(C);

    geo.connect(A, LocationType::Left, B);
    geo.connect(A, LocationType::Up, C);
    /*
        geo.connect(A, LocationType::Left, B, S_AB);    //Add Schur matrix here
    */

    A.set_boundary(LocationType::Right, PDEBoundaryType::Dirichlet);
    A.set_boundary(LocationType::Down, PDEBoundaryType::Dirichlet);

    B.set_boundary(LocationType::Down, PDEBoundaryType::Dirichlet);
    B.set_boundary(LocationType::Up, PDEBoundaryType::Dirichlet);
    B.set_boundary(LocationType::Left, PDEBoundaryType::Dirichlet);

    C.set_boundary(LocationType::Up, PDEBoundaryType::Dirichlet);
    C.set_boundary(LocationType::Left, PDEBoundaryType::Dirichlet);
    C.set_boundary(LocationType::Right, PDEBoundaryType::Dirichlet);

    geo.check();
    
    Variable p("p");
    field2 p_A;
    field2 p_B("p_B");
    field2 p_C("p_C");

    p.set_geometry(geo);
    p.set_field(A, p_A);
    p.set_field(B, p_B);
    p.set_field(C, p_C);

    std::cout << "main domain of p: " << p.geometry->main_domain->name << std::endl;
    std::cout <<  "name of p_A: " << p_A.get_name() << std::endl;

    // ConcatSolver2D solver(p);
    // solver.solve();

    return 0;
}
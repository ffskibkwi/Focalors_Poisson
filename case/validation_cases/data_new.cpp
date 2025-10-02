#include "core/domain/geometry2d.h"
#include "core/domain/domain2d.h"
#include "core/domain/variable.h"
#include "core/base/field2.h"

#include "core/boundary/boundary_type.h"

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

    geo.add_domain(A);
    geo.add_domain(B);
    geo.add_domain(C);

    geo.connect(A, LocationType::Left, B);
    geo.connect(A, LocationType::Up, C);

    A.set_boundary(LocationType::Up, PDEBoundaryType::Dirichlet);
    A.set_boundary(LocationType::Down, PDEBoundaryType::Dirichlet);

    geo.check();
    
    Variable p;
    field2 p_A("p_A");
    field2 p_B("p_B");
    field2 p_C("p_C");

    p.set_geometry(geo);
    p.set_field(A, p_A);    //In this step, must check if A belongs to geo
    p.set_field(B, p_B);
    p.set_field(C, p_C);
    
    // geo.solve_prepare();

    // void geometry::solve_prepare()
    // {
    //     //Find main domain
    //     //Find the domain has more than one neibour, that is the main domain
    //     //If there are more than one domain, the current version does not support
    // }
}
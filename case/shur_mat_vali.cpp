#include "pe/concat/domain_solver.h"
#include "pe/concat/schur_mat2d.h"
#include <iostream>
#include <string>

class DomainSolver2DTest : public DomainSolver2D
{
public:
    DomainSolver2DTest() {}
    ~DomainSolver2DTest() {}
    void solve(field2& f, bool is_debugmode = true)
    {
        for (int i = 0; i < f.get_nx(); i++)
        {
            for (int j = 0; j < f.get_ny(); j++)
            {
                f(i, j) *= i * f.get_ny() + j + 1;
            }
        }
    }
    bool is_comm_root() const { return true; }

    double get_hx() const { return 1; };
    double get_hy() const { return 1; };
};

template<typename SchurType>
void test_schur_direction(const std::string& label, Domain2DUniform& neighbor_domain, int nx, int ny)
{
    std::cout << "\n========== Testing " << label << " ==========" << std::endl;

    field2 f(nx, ny);
    for (int i = 0; i < f.get_nx(); i++)
    {
        for (int j = 0; j < f.get_ny(); j++)
        {
            f(i, j) = i * f.get_ny() + j + 1;
        }
    }

    std::cout << "Original Field:" << std::endl;
    f.print();

    SchurType          shur(&neighbor_domain);
    DomainSolver2DTest solver;
    shur.construct(&solver);

    std::cout << "Schur Matrix Structure:" << std::endl;
    shur.print();

    field2 result = shur * f;

    std::cout << "Result Field (after Schur * f):" << std::endl;
    result.print();
}

int main()
{
    int    nx          = 4;
    int    ny          = 3;
    int    neighbor_nx = 6;
    int    neighbor_ny = 3;
    double lx          = 1.0;
    double ly          = 1.0;

    Domain2DUniform neighbor_domain_x(neighbor_nx, neighbor_ny, lx, ly, "TestDomain");
    Domain2DUniform neighbor_domain_y(neighbor_ny, neighbor_nx, ly, lx, "TestDomain");

    test_schur_direction<SchurMat2D_left>("LEFT", neighbor_domain_x, nx, ny);
    test_schur_direction<SchurMat2D_right>("RIGHT", neighbor_domain_x, nx, ny);
    test_schur_direction<SchurMat2D_up>("UP", neighbor_domain_y, ny, nx);
    test_schur_direction<SchurMat2D_down>("DOWN", neighbor_domain_y, ny, nx);

    return 0;
}
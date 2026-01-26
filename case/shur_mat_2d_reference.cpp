#include "pe/concat/domain_solver.h"
#include "pe/concat/schur_mat2d.h"
#include <iostream>
#include <string>

class DomainSolver2DTest : public DomainSolver2D
{
public:
    DomainSolver2DTest() {}
    ~DomainSolver2DTest() {}
    void solve(field2& f)
    {
        for (int i = 0; i < f.get_nx(); i++)
        {
            for (int j = 0; j < f.get_ny(); j++)
            {
                f(i, j) *= i * f.get_ny() + j + 1;
            }
        }
    }
};

template<typename SchurType>
void test(const std::string& label, Domain2DUniform& neighbor_domain, int nx, int ny)
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
    int nx          = 4;
    int ny          = 3;
    int neighbor_nx = 6;
    int neighbor_ny = 3;

    Domain2DUniform neighbor_domain_x(neighbor_nx, neighbor_ny, "neighbor_domain_x");
    Domain2DUniform neighbor_domain_y(neighbor_ny, neighbor_nx, "neighbor_domain_y");

    test<SchurMat2D_left>("LEFT", neighbor_domain_x, nx, ny);
    test<SchurMat2D_right>("RIGHT", neighbor_domain_x, nx, ny);
    test<SchurMat2D_up>("UP", neighbor_domain_y, ny, nx);
    test<SchurMat2D_down>("DOWN", neighbor_domain_y, ny, nx);

    return 0;
}
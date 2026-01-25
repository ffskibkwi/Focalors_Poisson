#include "pe/concat/domain_solver.h"
#include "pe/concat/schur_mat2d.h"

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

int main()
{
    int    nx          = 4;
    int    ny          = 3;
    int    neighbor_nx = 6;
    int    neighbor_ny = 3;
    double lx          = 1.0;
    double ly          = 1.0;

    Domain2DUniform neighbor_domain(neighbor_nx, neighbor_ny, lx, ly, "Test");

    field2 f(nx, ny);

    for (int i = 0; i < f.get_nx(); i++)
    {
        for (int j = 0; j < f.get_ny(); j++)
        {
            f(i, j) = i * f.get_ny() + j + 1;
        }
    }

    f.print();

    std::cout << "---------------" << std::endl;

    SchurMat2D_left    shur(&neighbor_domain);
    DomainSolver2DTest solver;
    shur.construct(&solver);
    f = shur * f;

    shur.print();

    std::cout << "---------------" << std::endl;

    f.print();
}
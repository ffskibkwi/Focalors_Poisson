#include "base/domain/domain2d.h"
#include "base/domain/geometry2d.h"
#include "base/domain/geometry_tree.hpp"
#include "base/domain/variable2d.h"
#include "base/field/field2.h"
#include "pe/concat/gmres_solver2d.h"

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

void test(const std::string& label, LocationType neighbor_loc)
{
    std::cout << "\n========== Testing " << label << " ==========" << std::endl;

    int nx          = 3;
    int ny          = 5;
    int neighbor_n  = 4;
    int neighbor_nx = nx, neighbor_ny = ny;
    if (neighbor_loc == LocationType::Left || neighbor_loc == LocationType::Right)
        neighbor_nx = neighbor_n;
    else if (neighbor_loc == LocationType::Down || neighbor_loc == LocationType::Up)
        neighbor_ny = neighbor_n;

    int    m       = 2;
    double tol     = 1e-3;
    int    maxIter = 2;

    Geometry2D geo_tee;

    Domain2DUniform root(nx, ny, "root");
    Domain2DUniform neighbor(neighbor_nx, neighbor_ny, "neighbor");
    field2          p_root(nx, ny, "p_root");

    GMRESSolver2D       gmres(&root, m, tol, maxIter);
    DomainSolver2DTest* test_solver = new DomainSolver2DTest();
    gmres.set_solver(test_solver);

    std::unordered_map<LocationType, Domain2DUniform*> adj;
    adj[neighbor_loc] = &neighbor;

    std::unordered_map<Domain2DUniform*, DomainSolver2D*> solver_map;
    solver_map[&root]     = &gmres;
    solver_map[&neighbor] = test_solver;

    gmres.schur_mat_construct(adj, solver_map);

    for (int i = 0; i < p_root.get_nx(); i++)
    {
        for (int j = 0; j < p_root.get_ny(); j++)
        {
            p_root(i, j) = i * p_root.get_ny() + j + 1;
        }
    }

    gmres.solve(p_root);

    p_root.print();
}

int main()
{
    test("LEFT", LocationType::Left);
    test("RIGHT", LocationType::Right);
    test("DOWN", LocationType::Down);
    test("UP", LocationType::Up);
    return 0;
}
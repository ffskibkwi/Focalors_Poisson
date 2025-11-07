#include <mpi.h>
#include <iostream>
#include <cmath>

#include "io/config.h"
#include "base/domain/domain2d.h"
#include "base/domain/variable.h"
#include "base/location_boundary.h"
#include "pe_mpi/poisson/mpi_poisson_solver2d.h"

int main(int argc, char** argv)
{
    MPI_Init(&argc, &argv);

    int world_rank = 0, world_size = 1;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    // 控制并行规模与起始世界进程（argv[1]=num_proc, argv[2]=start_rank）
    int num_proc = world_size;
    int start_rank = 0;
    if (argc >= 2)
    {
        num_proc = std::atoi(argv[1]);
        if (num_proc <= 0 || num_proc > world_size) num_proc = world_size;
    }
    if (argc >= 3)
    {
        start_rank = std::atoi(argv[2]);
        if (start_rank < 0) start_rank = 0;
        if (start_rank >= world_size) start_rank = 0;
        if (start_rank + num_proc > world_size) num_proc = world_size - start_rank;
        if (num_proc < 0) num_proc = 0;
    }

    // 构造计算域与变量（简单方形区域，Dirichlet=0 边界）
    const int    nx = 16;
    const int    ny = 12;
    const double lx = 1.0;
    const double ly = 1.0;

    Domain2DUniform domain(nx, ny, lx, ly, "demo_domain");

    Variable var("phi");
    var.set_boundary_type(&domain, LocationType::Left,  PDEBoundaryType::Dirichlet);
    var.set_boundary_type(&domain, LocationType::Right, PDEBoundaryType::Dirichlet);
    var.set_boundary_type(&domain, LocationType::Down,  PDEBoundaryType::Dirichlet);
    var.set_boundary_type(&domain, LocationType::Up,    PDEBoundaryType::Dirichlet);
    var.set_boundary_value(&domain, LocationType::Left,  0.0);
    var.set_boundary_value(&domain, LocationType::Right, 0.0);
    var.set_boundary_value(&domain, LocationType::Down,  0.0);
    var.set_boundary_value(&domain, LocationType::Up,    0.0);

    // 环境配置（仅子通信器根打印求解阶段）
    EnvironmentConfig env;
    env.showCurrentStep = (world_rank == start_rank);

    // 仅子通信器根（world_rank==start_rank）持有全局 field
    field2 f_global;
    if (world_rank == start_rank)
    {
        domain.construct_field(f_global);
        // 构造一个简单的右端项 f(x,y)=1（常数源项），以测试端到端流程
        for (int i = 0; i < nx; ++i)
        {
            for (int j = 0; j < ny; ++j)
            {
                f_global(i, j) = 1.0;
            }
        }
    }

    // 构造 MPI 版本 Poisson 求解器（从 Variable 读取边界类型, 指定起始世界进程）
    MPIPoissonSolver2D solver(&domain, &var, num_proc, start_rank, &env, MPI_COMM_WORLD);

    // 求解：内部会在子通信器前 num_proc 个进程中并行，且自动完成分发/回收
    solver.solve(f_global);

    // 打印结果：仅子通信器根检查部分数据
    if (world_rank == start_rank)
    {
        std::cout << "[mpi_demo] Solution sample (5 points):\n";
        std::cout << "phi(0,0)=" << f_global(0,0)
                  << ", phi(nx/2,ny/2)=" << f_global(nx/2, ny/2)
                  << ", phi(nx-1,ny-1)=" << f_global(nx-1, ny-1) << "\n";
        std::cout << "phi(1,1)=" << f_global(1,1)
                  << ", phi(nx-2,ny-2)=" << f_global(nx-2, ny-2) << "\n";
    }

    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();
    return 0;
}



#pragma once

#include "base/pch.h"

#include "base/domain/domain2d.h"
#include "base/location_boundary.h"
#include "io/config.h"
#include "mpi_poisson_solver2d.h"
#include "pe/concat/schur_mat2d.h"
#include "pe/concat/domain_solver.h"
#include "pe/concat/gmres.h"


#include <mpi.h>
#include <unordered_map>

/**
 * @brief MPI 版本 GMRES 求解器（二维），内部 Poisson 使用 MPIPoissonSolver2D。
 *
 * 约定：右端项 b 仅由传入通信器的 rank 0 持有；内部 Poisson 求解器沿用同样约定。
 */
class MPIGMRESSolver2D : public DomainSolver2D
{
public:
    MPIGMRESSolver2D(Domain2DUniform*   in_domain,
                     Variable*          in_variable,
                     int                in_m,
                     double             in_tol,
                     int                in_maxIter,
                     MPI_Comm           in_comm,
                     EnvironmentConfig* in_env_config = nullptr);
    ~MPIGMRESSolver2D();

    void solve(field2& f, bool is_debugmode = true) override;
    // void solve_collective_root_owned(field2& f, bool is_debugmode = true) override;
    bool is_comm_root() const override { return comm_rank == 0; }

    double get_hx() const override { return domain->hx; }
    double get_hy() const override { return domain->hy; }

    void schur_mat_construct(const std::unordered_map<LocationType, Domain2DUniform*>&    adjacency_key,
                             const std::unordered_map<Domain2DUniform*, DomainSolver2D*>& solver_map);

    void set_initial_guess(field2* x0) { x0_override = x0; }

private:
    // 输入与配置
    Domain2DUniform*   domain;
    Variable*          variable   = nullptr;
    EnvironmentConfig* env_config = nullptr;
    EnvironmentConfig  inner_env_config; // 内部 Poisson：关闭 step 打印

    // Krylov 配置
    int    m       = 0;
    double tol     = 0.0;
    int    maxIter = 0;

    // 通信器
    MPI_Comm comm      = MPI_COMM_NULL;
    int      comm_rank = -1, comm_size = 0;

    // Schur 矩阵参数
    std::vector<SchurMat2D*> S_params;

    // 内部 Poisson（MPI）
    MPIPoissonSolver2D* pe_solver = nullptr;

    // 预分配缓冲
    std::vector<field2> V;
    std::vector<double> H, cs, sn, g, y, resVec;
    field2              x_buf, r_buf, w_buf, ft_buf, mul_buf, afun_buf;

    // 可选初值覆盖
    field2* x0_override = nullptr;

    // 内部方法
    field2& Afun(field2& x);
    void    maybe_print_res() const;
};

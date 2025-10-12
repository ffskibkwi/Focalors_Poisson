#pragma once

#include "pch.h"

#include "core/boundary/boundary_type.h"
#include "core/domain/domain2d.h"
#include "pe/poisson_new/poisson_solver2d.h"
#include "chasing_method2d.h"
#include "domain_solver.h"
#include "gmres.h"
#include "Schur_mat.h"
#include "pe/poisson_new/poisson_solver2d.h"
#include <unordered_map>

class GMRESSolver2D : public DomainSolver2D
{
public:
    GMRESSolver2D(Domain2DUniform*              in_domain,
                  int                           in_m,
                  double                        in_tol,
                  int                           in_maxIter);
    ~GMRESSolver2D();

    void solve(field2& f) override;

    void schur_mat_construct(const std::unordered_map<LocationType, Domain2DUniform*>& adjacency_key, 
                             const std::unordered_map<Domain2DUniform*, DomainSolver2D*>& solver_map);

    // 可选：设置外部初始猜测 x0；若未设置则默认 x0 = b
    void set_initial_guess(field2* x0) { x0_override = x0; }

private:
    // 成员参数（原 gmres 的未指定入参）：
    Domain2DUniform&              domain;
    std::vector<Schur_mat*>       S_params;
    int                           m = 0;
    double                        tol = 0.0;
    int                           maxIter = 0;

    PoissonSolver2D*              pe_solver;

    // 预分配的临时缓冲（避免在 solve 中动态分配）：
    std::vector<field2>           V;              // 大小 m+1，每个与域尺寸一致
    std::vector<double>           H;              // 大小 (m+1)*m
    std::vector<double>           cs;             // 大小 m
    std::vector<double>           sn;             // 大小 m
    std::vector<double>           g;              // 大小 m+1
    std::vector<double>           y;              // 最大大小 m（使用前 j 元素）
    std::vector<double>           resVec;         // 残差历史（每次 solve 前清空）

    // 复用的 field2 缓冲，避免 solve 内部动态分配
    field2                        x_buf;
    field2                        r_buf;
    field2                        w_buf;
    field2                        ft_buf;
    field2                        mul_buf;
    field2                        afun_buf;

    // 可选初值覆盖
    field2*                       x0_override = nullptr;

    // 与原实现一致的 Afun
    field2& Afun(field2& x);

    // field2 gmres(field2&                   b,
    //     field2&                   x,
    //     std::vector<Schur_mat*>&   S_params,
    //     PoissonSolver2DInterface& pe_solver,
    //     int                       m,
    //     double                    tol,
    //     int                       maxIter,
    //     std::vector<double>&      resVec);
};
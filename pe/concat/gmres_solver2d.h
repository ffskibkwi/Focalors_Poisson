#pragma once

#include "base/pch.h"

#include "base/config.h"
#include "base/domain/domain2d.h"
#include "base/location_boundary.h"
#include "domain_solver.h"
#include "schur_mat2d.h"
#include <unordered_map>

class GMRESSolver2D : public DomainSolver2D
{
public:
    GMRESSolver2D(Domain2DUniform* in_domain, int in_m, double in_tol, int in_maxIter);
    ~GMRESSolver2D();

    void set_solver(DomainSolver2D* _solver) { pe_solver = _solver; }
    void solve(field2& f) override;

    void schur_mat_construct(const std::unordered_map<LocationType, Domain2DUniform*>&    adjacency_key,
                             const std::unordered_map<Domain2DUniform*, DomainSolver2D*>& solver_map);

    // 可选：设置外部初始猜测 x0；若未设置则默认 x0 = b
    void set_initial_guess(field2* x0) { x0_override = x0; }

private:
    Domain2DUniform*         domain;
    std::vector<SchurMat2D*> S_params;
    int                      m       = 0;
    double                   tol     = 0.0;
    int                      maxIter = 0;

    DomainSolver2D* pe_solver;

    std::vector<field2> V;      // 大小 m+1，每个与域尺寸一致
    std::vector<double> H;      // 大小 (m+1)*m
    std::vector<double> cs;     // 大小 m
    std::vector<double> sn;     // 大小 m
    std::vector<double> g;      // 大小 m+1
    std::vector<double> y;      // 最大大小 m（使用前 j 元素）
    std::vector<double> resVec; // 残差历史（每次 solve 前清空）

    field2 x_buf;
    field2 r_buf;
    field2 w_buf;
    field2 ft_buf;
    field2 mul_buf;
    field2 afun_buf;

    field2* x0_override = nullptr;

    field2& Afun(field2& x);

    void maybe_print_res() const;

    int solve_call_count = 0;
};
#pragma once

#include "base/pch.h"

#include "base/domain/domain3d.h"
#include "base/location_boundary.h"
#include "domain_solver.h"
#include "io/config.h"
#include "schur_mat3d.h"
#include <unordered_map>

class GMRESSolver3D : public DomainSolver3D
{
public:
    GMRESSolver3D(Domain3DUniform*   in_domain,
                  int                in_m,
                  double             in_tol,
                  int                in_maxIter,
                  EnvironmentConfig* in_env_config = nullptr);
    ~GMRESSolver3D();

    void set_solver(DomainSolver3D* _solver) { pe_solver = _solver; }
    void solve(field3& f) override;

    void schur_mat_construct(const std::unordered_map<LocationType, Domain3DUniform*>&    adjacency_key,
                             const std::unordered_map<Domain3DUniform*, DomainSolver3D*>& solver_map);

    // 可选：设置外部初始猜测 x0；若未设置则默认 x0 = b
    void set_initial_guess(field3* x0) { x0_override = x0; }

private:
    // 成员参数
    Domain3DUniform*         domain;
    std::vector<SchurMat3D*> S_params;
    int                      m       = 0;
    double                   tol     = 0.0;
    int                      maxIter = 0;

    DomainSolver3D* pe_solver;

    std::vector<field3> V;      // 大小 m+1，每个与域尺寸一致
    std::vector<double> H;      // 大小 (m+1)*m
    std::vector<double> cs;     // 大小 m
    std::vector<double> sn;     // 大小 m
    std::vector<double> g;      // 大小 m+1
    std::vector<double> y;      // 最大大小 m（使用前 j 元素）
    std::vector<double> resVec; // 残差历史（每次 solve 前清空）

    field3 x_buf;
    field3 r_buf;
    field3 w_buf;
    field3 ft_buf;
    field3 mul_buf;
    field3 afun_buf;

    field3* x0_override = nullptr;

    EnvironmentConfig* env_config = nullptr;

    field3& Afun(field3& x);

    void maybe_print_res() const;
};
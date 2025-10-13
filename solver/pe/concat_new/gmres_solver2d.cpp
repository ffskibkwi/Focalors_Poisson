#include "gmres_solver2d.h"
#include <cmath>
#include <vector>

GMRESSolver2D::GMRESSolver2D(Domain2DUniform*              in_domain,
                             int                           in_m,
                             double                        in_tol,
                             int                           in_maxIter)
    : domain(in_domain)
    , m(in_m)
    , tol(in_tol)
    , maxIter(in_maxIter)
{
    // 预分配 field2 缓冲与 Krylov 基
    int nx = domain->nx;
    int ny = domain->ny;

    pe_solver = new PoissonSolver2D(domain);

    x_buf.init(nx, ny);
    r_buf.init(nx, ny);
    w_buf.init(nx, ny);
    ft_buf.init(nx, ny);
    mul_buf.init(nx, ny);
    afun_buf.init(nx, ny);

    V.resize(m + 1);
    for (int i = 0; i <= m; i++)
        V[i].init(domain->nx, domain->ny);

    H.resize((m + 1) * m, 0.0);
    cs.resize(m, 0.0);
    sn.resize(m, 0.0);
    g.resize(m + 1, 0.0);
    y.resize(m, 0.0);
    resVec.clear();
}

GMRESSolver2D::~GMRESSolver2D()
{
    delete pe_solver;
}

void GMRESSolver2D::schur_mat_construct(const std::unordered_map<LocationType, Domain2DUniform*>& adjacency_key, const std::unordered_map<Domain2DUniform*, DomainSolver2D*>& solver_map)
{
    for (auto &[location, neighbour_domain] : adjacency_key)
    {
        //Construct the Schur matrix for each neibour domain of main domain
        Schur_mat* current = nullptr;
        switch (location)
        {
            case LocationType::Left:
            {
                current = new Schur_mat_left(*domain, *neighbour_domain);
                current->construct(solver_map.at(neighbour_domain));
                S_params.push_back(current);
            }  
            break;
            case LocationType::Right:
            {
                current = new Schur_mat_right(*domain, *neighbour_domain);
                current->construct(solver_map.at(neighbour_domain));
                S_params.push_back(current);
            }
            break;
            case LocationType::Up:
            {
                current = new Schur_mat_up(*domain, *neighbour_domain);
                current->construct(solver_map.at(neighbour_domain));
                S_params.push_back(current);
            }
            break;
            case LocationType::Down:
            {
                current = new Schur_mat_down(*domain, *neighbour_domain);
                current->construct(solver_map.at(neighbour_domain));
                S_params.push_back(current);
            }
            break;
            default:
                throw std::invalid_argument("Invalid location type");
        }
    }
}
field2& GMRESSolver2D::Afun(field2& x)
{
    // 与 concat_new/gmres.cpp 一致：ft = sum(S_i * x); pe_solver->solve(ft); return x - ft
    ft_buf = x;
    ft_buf.clear(0.);
    for (auto& s : S_params)
    {
        mul_buf = (*s) * x;                // 仍会产生一次临时，但接口所限
        ft_buf.add_affine_transform(1.0, mul_buf, 0.0); // ft += mul_buf
    }
    pe_solver->solve(ft_buf);

    afun_buf = x;                          // afun_buf = x - ft_buf（避免新分配）
    afun_buf.add_affine_transform(-1.0, ft_buf, 0.0);
    return afun_buf;
}

void GMRESSolver2D::solve(field2& b)
{
    // 每次执行前清空残差历史
    resVec.clear();

    // 初始解 x0：若未设置则默认 x0 = b
    if (x0_override)
        x_buf = *x0_override;
    else
        x_buf = b;

    double beta = 0.0;

    for (int outer = 0; outer < maxIter; outer++)
    {
        {
            // r = b - Afun(x)
            field2& Ax = Afun(x_buf);
            r_buf = b;
            r_buf.add_affine_transform(-1.0, Ax, 0.0);
            beta = r_buf.norm();
            resVec.push_back(beta);
            if (beta < tol)
            {
                b = x_buf; // 直接覆盖 b
                return;
            }
        }

        // 重置 H, cs, sn, g
        std::fill(H.begin(), H.end(), 0.0);
        std::fill(cs.begin(), cs.end(), 0.0);
        std::fill(sn.begin(), sn.end(), 0.0);
        std::fill(g.begin(), g.end(), 0.0);

        // V[0] = r / beta
        V[0].clear(0.0);
        V[0].add_affine_transform(1.0 / beta, r_buf, 0.0);
        g[0] = beta;

        int j = 0;
        for (; j < m; j++)
        {
            // w = Afun(V[j])
            {
                field2& AVj = Afun(V[j]);
                w_buf = AVj;

                // Arnoldi 正交化
                for (int i = 0; i <= j; i++)
                {
                    H[i * m + j] = V[i].dot(w_buf);
                    w_buf.add_affine_transform(-H[i * m + j], V[i], 0.0); // w -= H(i,j)*V[i]
                }

                const double h_j1j = w_buf.norm();
                if (h_j1j < 1e-12)
                    break; // 提前终止

                H[(j + 1) * m + j] = h_j1j;
                V[j + 1].clear(0.0);
                V[j + 1].add_affine_transform(1.0 / h_j1j, w_buf, 0.0);
            }

            // 应用已有的 Givens 旋转到 H 的第 j 列
            for (int i = 0; i < j; i++)
            {
                const double temp       = cs[i] * H[i * m + j] + sn[i] * H[(i + 1) * m + j];
                H[(i + 1) * m + j]      = -sn[i] * H[i * m + j] + cs[i] * H[(i + 1) * m + j];
                H[i * m + j]            = temp;
            }

            // 新的 Givens 旋转
            const double h_jj  = H[j * m + j];
            const double h_j1j = H[(j + 1) * m + j];
            const double denom = std::hypot(h_jj, h_j1j);
            cs[j]              = h_jj / denom;
            sn[j]              = h_j1j / denom;

            // 更新 H 和 g
            H[j * m + j]       = cs[j] * H[j * m + j] + sn[j] * H[(j + 1) * m + j];
            H[(j + 1) * m + j] = 0.0;

            const double temp = cs[j] * g[j];
            g[j + 1]          = -sn[j] * g[j];
            g[j]              = temp;

            if (std::abs(g[j + 1]) < tol)
            {
                j++; // 包含当前 j
                break;
            }
        }

        // 回代求解 y（大小 j）
        for (int i = 0; i < m; i++)
            y[i] = 0.0;
        for (int i = j - 1; i >= 0; i--)
        {
            y[i] = g[i];
            for (int k = i + 1; k < j; k++)
                y[i] -= H[i * m + k] * y[k];
            y[i] /= H[i * m + i];
        }

        // 更新解 x = x + V(:,0:j-1)*y
        for (int i = 0; i < j; i++)
            x_buf.add_affine_transform(y[i], V[i], 0.0);
    }

    // 达到最大迭代，返回当前近似解
    b = x_buf;
}
#include "gmres_solver3d.h"

namespace
{
    const char* location_to_string(LocationType location)
    {
        switch (location)
        {
            case LocationType::Left:
                return "Left";
            case LocationType::Right:
                return "Right";
            case LocationType::Front:
                return "Front";
            case LocationType::Back:
                return "Back";
            case LocationType::Down:
                return "Down";
            case LocationType::Up:
                return "Up";
            default:
                return "Unknown";
        }
    }

    const char* safe_domain_name(const Domain3DUniform* domain) { return domain ? domain->name.c_str() : "unknown"; }

    void
    print_gmres_done(const EnvironmentConfig* env, const Domain3DUniform* domain, const std::vector<double>& res_vec)
    {
        if (!env || !env->showCurrentStep)
            return;

        std::cout << "[GMRES3D] solve: done (domain " << safe_domain_name(domain) << ")";
        if (!res_vec.empty())
            std::cout << " iter=" << res_vec.size() << " final_res=" << res_vec.back();
        else
            std::cout << " iter=0 final_res=n/a";
        std::cout << std::endl;
    }
} // namespace

GMRESSolver3D::GMRESSolver3D(Domain3DUniform*   in_domain,
                             Variable3D*        in_variable,
                             int                in_m,
                             double             in_tol,
                             int                in_maxIter,
                             EnvironmentConfig* in_env_config)
    : domain(in_domain)
    , variable(in_variable)
    , m(in_m)
    , tol(in_tol)
    , maxIter(in_maxIter)
{
    env_config = in_env_config;

    // 预分配 field3 缓冲与 Krylov 基
    int nx = domain->nx;
    int ny = domain->ny;
    int nz = domain->nz;

    x_buf.init(nx, ny, nz, "x_buf");
    r_buf.init(nx, ny, nz, "r_buf");
    w_buf.init(nx, ny, nz, "w_buf");
    ft_buf.init(nx, ny, nz, "ft_buf");
    mul_buf.init(nx, ny, nz, "mul_buf");
    afun_buf.init(nx, ny, nz, "afun_buf");

    V.resize(m + 1);
    for (int i = 0; i <= m; i++)
        V[i].init(nx, ny, nz, "V_" + std::to_string(i));

    H.resize((m + 1) * m, 0.0);
    cs.resize(m, 0.0);
    sn.resize(m, 0.0);
    g.resize(m + 1, 0.0);
    y.resize(m, 0.0);
    resVec.clear();
}

GMRESSolver3D::~GMRESSolver3D()
{
    delete pe_solver;
    for (auto* S : S_params)
        delete S;
}

void GMRESSolver3D::schur_mat_construct(const std::unordered_map<LocationType, Domain3DUniform*>&    adjacency_key,
                                        const std::unordered_map<Domain3DUniform*, DomainSolver3D*>& solver_map)
{
    if (env_config && env_config->showCurrentStep)
        std::cout << "[GMRES3D] Schur construct: domain " << safe_domain_name(domain)
                  << " neighbors=" << adjacency_key.size() << std::endl;
    for (auto& [location, neighbour_domain] : adjacency_key)
    {
        if (env_config && env_config->showCurrentStep)
        {
            std::cout << "[GMRES3D] Schur construct: " << safe_domain_name(domain) << " <-> "
                      << safe_domain_name(neighbour_domain) << " at " << location_to_string(location) << " start"
                      << std::endl;
        }
        // Construct the Schur matrix for each neighbour domain of main domain
        SchurMat3D* current = nullptr;
        switch (location)
        {
            case LocationType::Left: {
                current = new SchurMat3D_left(neighbour_domain);
                current->construct(solver_map.at(neighbour_domain));
                S_params.push_back(current);
            }
            break;
            case LocationType::Right: {
                current = new SchurMat3D_right(neighbour_domain);
                current->construct(solver_map.at(neighbour_domain));
                S_params.push_back(current);
            }
            break;
            case LocationType::Front: {
                current = new SchurMat3D_front(neighbour_domain);
                current->construct(solver_map.at(neighbour_domain));
                S_params.push_back(current);
            }
            break;
            case LocationType::Back: {
                current = new SchurMat3D_back(neighbour_domain);
                current->construct(solver_map.at(neighbour_domain));
                S_params.push_back(current);
            }
            break;
            case LocationType::Down: {
                current = new SchurMat3D_down(neighbour_domain);
                current->construct(solver_map.at(neighbour_domain));
                S_params.push_back(current);
            }
            break;
            case LocationType::Up: {
                current = new SchurMat3D_up(neighbour_domain);
                current->construct(solver_map.at(neighbour_domain));
                S_params.push_back(current);
            }
            break;
            default:
                throw std::invalid_argument("Invalid location type");
        }
        if (env_config && env_config->showCurrentStep)
        {
            std::cout << "[GMRES3D] Schur construct: " << safe_domain_name(domain) << " <-> "
                      << safe_domain_name(neighbour_domain) << " at " << location_to_string(location) << " done"
                      << std::endl;
        }
    }
    if (env_config && env_config->showCurrentStep)
        std::cout << "[GMRES3D] Schur construct: domain " << safe_domain_name(domain) << " done" << std::endl;
}

field3& GMRESSolver3D::Afun(field3& x)
{
    // ft = sum(S_i * x); pe_solver->solve(ft); return x - ft
    ft_buf.clear(0.);
    for (auto& s : S_params)
    {
        mul_buf = (*s) * x;                             // Schur 矩阵乘法
        ft_buf.add_affine_transform(1.0, mul_buf, 0.0); // ft += mul_buf
    }
    pe_solver->solve(ft_buf);

    afun_buf = x; // afun_buf = x - ft_buf
    afun_buf.add_affine_transform(-1.0, ft_buf, 0.0);
    return afun_buf;
}

void GMRESSolver3D::maybe_print_res() const
{
    if (env_config && env_config->showGmresRes)
    {
        std::cout << "GMRES3D resVec: [";
        for (size_t i = 0; i < resVec.size(); ++i)
        {
            std::cout << resVec[i];
            if (i + 1 < resVec.size())
                std::cout << ", ";
        }
        std::cout << "]" << std::endl;
    }
}

void GMRESSolver3D::solve(field3& b)
{
    if (env_config && env_config->showCurrentStep)
        std::cout << "[GMRES3D] solve: start (domain " << safe_domain_name(domain) << ")" << std::endl;

    // Actually the solver is for the equation (I-{A^-1}S)x={A^-1}b
    pe_solver->solve(b);

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
            field3& Ax = Afun(x_buf);
            r_buf      = b;
            r_buf.add_affine_transform(-1.0, Ax, 0.0);

            // 计算残差范数
            beta = 0.0;
            for (int i = 0; i < r_buf.get_nx(); i++)
                for (int j = 0; j < r_buf.get_ny(); j++)
                    for (int k = 0; k < r_buf.get_nz(); k++)
                        beta += r_buf(i, j, k) * r_buf(i, j, k);
            beta = std::sqrt(beta);

            resVec.push_back(beta);
            if (beta < tol)
            {
                b = x_buf; // 直接覆盖 b
                maybe_print_res();
                print_gmres_done(env_config, domain, resVec);
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
        double inv_beta = 1.0 / beta;
        for (int i = 0; i < r_buf.get_nx(); i++)
            for (int j = 0; j < r_buf.get_ny(); j++)
                for (int k = 0; k < r_buf.get_nz(); k++)
                    V[0](i, j, k) = r_buf(i, j, k) * inv_beta;
        g[0] = beta;

        int j = 0;
        for (; j < m; j++)
        {
            // w = Afun(V[j])
            {
                field3& AVj = Afun(V[j]);
                w_buf       = AVj;

                // Arnoldi 正交化
                for (int i = 0; i <= j; i++)
                {
                    // 计算点积
                    double dot_product = 0.0;
                    for (int ii = 0; ii < w_buf.get_nx(); ii++)
                        for (int jj = 0; jj < w_buf.get_ny(); jj++)
                            for (int kk = 0; kk < w_buf.get_nz(); kk++)
                                dot_product += V[i](ii, jj, kk) * w_buf(ii, jj, kk);

                    H[i * m + j] = dot_product;
                    // w -= H(i,j)*V[i]
                    w_buf.add_affine_transform(-H[i * m + j], V[i], 0.0);
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
                const double temp  = cs[i] * H[i * m + j] + sn[i] * H[(i + 1) * m + j];
                H[(i + 1) * m + j] = -sn[i] * H[i * m + j] + cs[i] * H[(i + 1) * m + j];
                H[i * m + j]       = temp;
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
    maybe_print_res();
    print_gmres_done(env_config, domain, resVec);
}

#include "mpi_gmres_solver2d.h"
#include <climits>
#include <cmath>
#include <iostream>
#include <vector>

MPIGMRESSolver2D::MPIGMRESSolver2D(Domain2DUniform*   in_domain,
                                   Variable*          in_variable,
                                   int                in_m,
                                   double             in_tol,
                                   int                in_maxIter,
                                   MPI_Comm           in_comm,
                                   EnvironmentConfig* in_env_config)
    : domain(in_domain)
    , variable(in_variable)
    , env_config(in_env_config)
    , m(in_m)
    , tol(in_tol)
    , maxIter(in_maxIter)
    , comm(in_comm)
{
    MPI_Comm_rank(comm, &comm_rank);
    MPI_Comm_size(comm, &comm_size);

    // 内部 Poisson 的配置：关闭 step 打印，GMRES 残差信息由外层控制
    inner_env_config = EnvironmentConfig {};
    if (env_config)
        inner_env_config.showGmresRes = env_config->showGmresRes;
    inner_env_config.showCurrentStep = false;

    // 内部 Poisson（MPI）遵循“仅 comm rank0 持有 b”的约定
    pe_solver = new MPIPoissonSolver2D(domain, variable, comm_size, /*start_rank*/ 0, &inner_env_config, comm);

    // 预分配 field2 缓冲与 Krylov 基
    int nx = domain->nx;
    int ny = domain->ny;
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

MPIGMRESSolver2D::~MPIGMRESSolver2D() { delete pe_solver; }

void MPIGMRESSolver2D::solve_collective_root_owned(field2& f, bool is_debugmode)
{
    if (comm_rank == 0)
        solve(f, is_debugmode);
    else
    {
        field2 dummy;
        dummy.init(domain->nx, domain->ny);
        solve(dummy, is_debugmode);
    }
}

void MPIGMRESSolver2D::schur_mat_construct(const std::unordered_map<LocationType, Domain2DUniform*>&    adjacency_key,
                                           const std::unordered_map<Domain2DUniform*, DomainSolver2D*>& solver_map)
{
    if (env_config && env_config->showCurrentStep && comm_rank == 0)
        std::cout << "[MPIGMRES] Schur construct: start" << std::endl;
    for (auto& [location, neighbour_domain] : adjacency_key)
    {
        Schur_mat* current = nullptr;
        switch (location)
        {
            case LocationType::Left:
                current = new Schur_mat_left(*domain, *neighbour_domain);
                break;
            case LocationType::Right:
                current = new Schur_mat_right(*domain, *neighbour_domain);
                break;
            case LocationType::Up:
                current = new Schur_mat_up(*domain, *neighbour_domain);
                break;
            case LocationType::Down:
                current = new Schur_mat_down(*domain, *neighbour_domain);
                break;
            default:
                throw std::invalid_argument("Invalid location type");
        }

        // 仅在拥有子域 solver 的进程上参与 Schur 构造（集体调用 child's solve），
        // 然后由“同时是子通信器 root 的父通信器进程”打包并在父通信器内广播。
        bool has_child_solver = false;
        {
            auto it          = solver_map.find(neighbour_domain);
            has_child_solver = (it != solver_map.end() && it->second != nullptr);
        }
        int root_candidate =
            (has_child_solver && solver_map.at(neighbour_domain)->is_comm_root()) ? comm_rank : INT_MAX;
        int root_rank = INT_MAX;
        MPI_Allreduce(&root_candidate, &root_rank, 1, MPI_INT, MPI_MIN, comm);
        if (root_rank == INT_MAX)
        {
            if (comm_rank == 0)
                std::cerr << "[MPIGMRES] Schur construct: no available child solver for location, domain="
                          << neighbour_domain->name << std::endl;
            delete current;
            throw std::runtime_error("MPIGMRES Schur construct failed: no child solver present in parent communicator");
        }

        int                 cosize = current->get_size();
        std::vector<double> matbuf(static_cast<size_t>(cosize) * static_cast<size_t>(cosize), 0.0);

        if (has_child_solver)
        {
            DomainSolver2D* child_solver = solver_map.at(neighbour_domain);
            current->construct(child_solver);
            if (comm_rank == root_rank)
            {
                // 打包为行主序
                current->store_rowmajor(matbuf.data());
            }
        }
        // 广播 Schur 矩阵到父通信器所有进程
        MPI_Bcast(matbuf.data(), static_cast<int>(matbuf.size()), MPI_DOUBLE, root_rank, comm);
        // 反填充到当前进程的 current 中
        if (comm_rank != root_rank)
            current->load_rowmajor(matbuf.data());

        S_params.push_back(current);
    }
    if (env_config && env_config->showCurrentStep && comm_rank == 0)
        std::cout << "[MPIGMRES] Schur construct: done" << std::endl;
}

field2& MPIGMRESSolver2D::Afun(field2& x)
{
    // 与串行一致：ft = sum(S_i * x); pe_solver->solve(ft); return x - ft
    ft_buf = x;
    ft_buf.clear(0.);
    for (auto& s : S_params)
    {
        mul_buf = (*s) * x;
        ft_buf.add_affine_transform(1.0, mul_buf, 0.0);
    }
    pe_solver->solve(ft_buf, false);

    afun_buf = x;
    afun_buf.add_affine_transform(-1.0, ft_buf, 0.0);
    return afun_buf;
}

void MPIGMRESSolver2D::maybe_print_res() const
{
    if (env_config && env_config->showGmresRes && comm_rank == 0)
    {
        std::cout << "GMRES resVec: [";
        for (size_t i = 0; i < resVec.size(); ++i)
        {
            std::cout << resVec[i];
            if (i + 1 < resVec.size())
                std::cout << ", ";
        }
        std::cout << "]" << std::endl;
    }
}

void MPIGMRESSolver2D::solve(field2& b, bool is_debugmode)
{
    if (env_config && env_config->showCurrentStep && comm_rank == 0)
        std::cout << "[MPIGMRES] solve: start" << std::endl;

    // 实际上求解 (I - A^{-1}S) x = A^{-1} b
    if (env_config && env_config->showCurrentStep && comm_rank == 0)
        std::cout << "[MPIGMRES] initial Poisson: start" << std::endl;
    pe_solver->solve(b, is_debugmode);
    if (env_config && env_config->showCurrentStep && comm_rank == 0)
        std::cout << "[MPIGMRES] initial Poisson: done" << std::endl;

    // 清空残差历史
    resVec.clear();

    // 初始解 x0：若未设置则默认 x0 = b
    if (x0_override)
        x_buf = *x0_override;
    else
        x_buf = b;

    double beta = 0.0;
    // 一致化 beta（防止积累误差导致分支不一致）
    beta = r_buf.norm(); // 将在第一轮内赋值，这里先置 0

    for (int outer = 0; outer < maxIter; outer++)
    {
        // r = b - Afun(x)
        {
            field2& Ax = Afun(x_buf);
            r_buf      = b;
            r_buf.add_affine_transform(-1.0, Ax, 0.0);
            beta = r_buf.norm();
            // 各进程一致化 beta（取最大或平均均可，这里取最大更稳定）
            double beta_global = 0.0;
            MPI_Allreduce(&beta, &beta_global, 1, MPI_DOUBLE, MPI_MAX, comm);
            beta = beta_global;
            resVec.push_back(beta);
            if (beta < tol)
            {
                b = x_buf;
                maybe_print_res();
                return;
            }
        }

        std::fill(H.begin(), H.end(), 0.0);
        std::fill(cs.begin(), cs.end(), 0.0);
        std::fill(sn.begin(), sn.end(), 0.0);
        std::fill(g.begin(), g.end(), 0.0);

        V[0].clear(0.0);
        V[0].add_affine_transform(1.0 / beta, r_buf, 0.0);
        g[0] = beta;

        int j = 0;
        for (; j < m; j++)
        {
            // w = Afun(V[j])
            {
                if (env_config && env_config->showCurrentStep && comm_rank == 0)
                    std::cout << "[MPIGMRES] Afun Poisson: start (j=" << j << ")" << std::endl;
                field2& AVj = Afun(V[j]);
                if (env_config && env_config->showCurrentStep && comm_rank == 0)
                    std::cout << "[MPIGMRES] Afun Poisson: done  (j=" << j << ")" << std::endl;
                w_buf = AVj;

                // Arnoldi
                for (int i = 0; i <= j; i++)
                {
                    H[i * m + j] = V[i].dot(w_buf);
                    w_buf.add_affine_transform(-H[i * m + j], V[i], 0.0);
                }

                double h_j1j = w_buf.norm();
                // 各进程一致化 h_j1j
                double h_global = 0.0;
                MPI_Allreduce(&h_j1j, &h_global, 1, MPI_DOUBLE, MPI_MAX, comm);
                h_j1j = h_global;
                if (h_j1j < 1e-12)
                    break;

                H[(j + 1) * m + j] = h_j1j;
                V[j + 1].clear(0.0);
                V[j + 1].add_affine_transform(1.0 / h_j1j, w_buf, 0.0);
            }

            // 应用已有的 Givens
            for (int i = 0; i < j; i++)
            {
                const double temp  = cs[i] * H[i * m + j] + sn[i] * H[(i + 1) * m + j];
                H[(i + 1) * m + j] = -sn[i] * H[i * m + j] + cs[i] * H[(i + 1) * m + j];
                H[i * m + j]       = temp;
            }

            // 新的 Givens
            const double h_jj  = H[j * m + j];
            const double h_j1j = H[(j + 1) * m + j];
            const double denom = std::hypot(h_jj, h_j1j);
            cs[j]              = h_jj / denom;
            sn[j]              = h_j1j / denom;

            H[j * m + j]       = cs[j] * H[j * m + j] + sn[j] * H[(j + 1) * m + j];
            H[(j + 1) * m + j] = 0.0;

            const double temp = cs[j] * g[j];
            g[j + 1]          = -sn[j] * g[j];
            g[j]              = temp;

            double gj1_abs        = std::abs(g[j + 1]);
            double gj1_abs_global = 0.0;
            MPI_Allreduce(&gj1_abs, &gj1_abs_global, 1, MPI_DOUBLE, MPI_MAX, comm);
            if (gj1_abs_global < tol)
            {
                j++;
                break;
            }
        }

        for (int i = 0; i < m; i++)
            y[i] = 0.0;
        for (int i = j - 1; i >= 0; i--)
        {
            y[i] = g[i];
            for (int k = i + 1; k < j; k++)
                y[i] -= H[i * m + k] * y[k];
            y[i] /= H[i * m + i];
        }

        for (int i = 0; i < j; i++)
            x_buf.add_affine_transform(y[i], V[i], 0.0);
    }

    b = x_buf;
    maybe_print_res();
    if (env_config && env_config->showCurrentStep && comm_rank == 0)
        std::cout << "[MPIGMRES] solve: done" << std::endl;
}

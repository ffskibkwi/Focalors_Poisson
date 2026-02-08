#include "gmres_solver2d_slab_x.h"
#include "base/parallel/mpi/mpi_misc.h"
#include "instrumentor/timer.h"
#include "io/csv_writer_2d.h"

GMRESSolver2DSlabX::GMRESSolver2DSlabX(Domain2DUniform* in_domain,
                                       int              in_m,
                                       double           in_tol,
                                       int              in_maxIter,
                                       MPI_Comm         _communicator)
    : domain(in_domain)
    , m(in_m)
    , tol(in_tol)
    , maxIter(in_maxIter)
    , communicator(_communicator)
{
    // 预分配 field2 缓冲与 Krylov 基
    int nx = domain->nx;
    int ny = domain->ny;

    MPI_Comm_rank(communicator, &mpi_rank);
    MPI_Comm_size(communicator, &mpi_size);

    int snx = MPIUtils::get_slab_length(nx, mpi_rank, mpi_size);

    x_buf.init(snx, ny);
    r_buf.init(snx, ny);
    w_buf.init(snx, ny);
    ft_buf.init(snx, ny);
    mul_buf.init(snx, ny);
    afun_buf.init(snx, ny);

    V.resize(m + 1);
    for (int i = 0; i <= m; i++)
        V[i].init(snx, ny);

    H.resize((m + 1) * m, 0.0);
    cs.resize(m, 0.0);
    sn.resize(m, 0.0);
    g.resize(m + 1, 0.0);
    y.resize(m, 0.0);
    resVec.clear();
}

GMRESSolver2DSlabX::~GMRESSolver2DSlabX()
{
    delete pe_solver;
    for (auto* S : S_params)
        delete S;
}

void GMRESSolver2DSlabX::schur_mat_construct(const std::unordered_map<LocationType, Domain2DUniform*>&    adjacency_key,
                                             const std::unordered_map<Domain2DUniform*, DomainSolver2D*>& solver_map)
{
    EnvironmentConfig& env_cfg = EnvironmentConfig::Get();

    if (env_cfg.showCurrentStep)
        std::cout << "[GMRES] Schur construct: start" << std::endl;
    for (auto& [location, neighbour_domain] : adjacency_key)
    {
        // Construct the Schur matrix for each neibour domain of main domain
        SchurMat2DSlabX* current = nullptr;

        DomainSolver2D* branch_solver = solver_map.at(neighbour_domain);

        switch (location)
        {
            case LocationType::Left: {
                current = new SchurMat2DSlabX_left(neighbour_domain, communicator);
                current->set_name("S_" + domain->name + "_Left_" + neighbour_domain->name);
                current->construct(branch_solver);
                if (env_cfg.debug_gmres)
                    current->write_csv(env_cfg.debugOutputDir);
                S_params.push_back(current);
            }
            break;
            case LocationType::Right: {
                current = new SchurMat2DSlabX_right(neighbour_domain, communicator);
                current->set_name("S_" + domain->name + "_Right_" + neighbour_domain->name);
                current->construct(branch_solver);
                if (env_cfg.debug_gmres)
                    current->write_csv(env_cfg.debugOutputDir);
                S_params.push_back(current);
            }
            break;
            case LocationType::Up: {
                current = new SchurMat2DSlabX_up(neighbour_domain, communicator);
                current->set_name("S_" + domain->name + "_Up_" + neighbour_domain->name);
                current->construct(branch_solver);
                if (env_cfg.debug_gmres)
                    current->write_csv(env_cfg.debugOutputDir);
                S_params.push_back(current);
            }
            break;
            case LocationType::Down: {
                current = new SchurMat2DSlabX_down(neighbour_domain, communicator);
                current->set_name("S_" + domain->name + "_Down_" + neighbour_domain->name);
                current->construct(branch_solver);
                if (env_cfg.debug_gmres)
                    current->write_csv(env_cfg.debugOutputDir);
                S_params.push_back(current);
            }
            break;
            default:
                throw std::invalid_argument("Invalid location type");
        }
    }
    if (env_cfg.showCurrentStep)
        std::cout << "[GMRES] Schur construct: done" << std::endl;
}
field2& GMRESSolver2DSlabX::Afun(field2& x)
{
    ft_buf = x;
    ft_buf.clear(0.);
    for (auto& s : S_params)
    {
        mul_buf = (*s) * x;                             // 仍会产生一次临时，但接口所限
        ft_buf.add_affine_transform(1.0, mul_buf, 0.0); // ft += mul_buf
    }
    pe_solver->solve(ft_buf);

    afun_buf = x; // afun_buf = x - ft_buf（避免新分配）
    afun_buf.add_affine_transform(-1.0, ft_buf, 0.0);
    return afun_buf;
}

void GMRESSolver2DSlabX::maybe_print_res() const
{
    EnvironmentConfig& env_cfg = EnvironmentConfig::Get();

    if (env_cfg.showGmresRes)
    {
        if (mpi_rank == 0)
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
}

void GMRESSolver2DSlabX::solve(field2& b)
{
    EnvironmentConfig& env_cfg = EnvironmentConfig::Get();

    SCOPE_TIMER("GMRESSolver2DSlabX::solve", TimeRecordType::None, env_cfg.track_pe_solve_detail_time);
    SCOPE_TIMER(env_cfg.pe_solve_total_name, TimeRecordType::Accumulate, false);

    if (env_cfg.showCurrentStep)
        std::cout << "[GMRES] solve: start" << std::endl;

    if (env_cfg.debug_gmres)
    {
        std::string fname_rhs =
            env_cfg.debugOutputDir + "/rhs_" + domain->name + "_" + std::to_string(solve_call_count);
        IO::write_csv(b, fname_rhs);
    }

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
            field2& Ax = Afun(x_buf);
            r_buf      = b;
            r_buf.add_affine_transform(-1.0, Ax, 0.0);

            beta                = r_buf.squared_sum();
            double beta_reduced = 0.0;
            MPI_Allreduce(&beta, &beta_reduced, 1, MPI_DOUBLE, MPI_SUM, communicator);
            beta = std::sqrt(beta_reduced);
            resVec.push_back(beta);

            if (beta < tol)
            {
                b = x_buf; // 直接覆盖 b
                maybe_print_res();
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
                w_buf       = AVj;

                // Arnoldi 正交化
                for (int i = 0; i <= j; i++)
                {
                    double H_temp = V[i].dot(w_buf);
                    MPI_Allreduce(&H_temp, &H[i * m + j], 1, MPI_DOUBLE, MPI_SUM, communicator);
                    w_buf.add_affine_transform(-H[i * m + j], V[i], 0.0); // w -= H(i,j)*V[i]
                }

                double h_temp    = w_buf.squared_sum();
                double h_reduced = 0.0;
                MPI_Allreduce(&h_temp, &h_reduced, 1, MPI_DOUBLE, MPI_SUM, communicator);
                const double h_j1j = std::sqrt(h_reduced);
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

    if (env_cfg.debug_gmres)
    {
        std::string fname_sol =
            env_cfg.debugOutputDir + "/sol_" + domain->name + "_" + std::to_string(solve_call_count);
        IO::write_csv(b, fname_sol);
    }
    solve_call_count++;

    if (env_cfg.showCurrentStep)
        std::cout << "[GMRES] solve: done" << std::endl;
}

void migrate_from(GMRESSolver2DSlabX* src, GMRESSolver2DSlabX* dest)
{
    if (dest != nullptr)
        dest->S_params.clear();

    int length_src = 0;
    if (src != nullptr)
        length_src = src->S_params.size();

    // sync to dest
    MPI_Allreduce(MPI_IN_PLACE, &length_src, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);

    // use int becuase the enum is uint8_t, unmatch with MPI_INT
    std::vector<int> loc_src(length_src, -1);
    if (src != nullptr)
    {
        for (int i = 0; i < length_src; i++)
        {
            SchurMat2DSlabX* S_src = src->S_params[i];
            loc_src[i]             = static_cast<int>(S_src->get_loc());
        }
    }

    // sync to dest
    MPI_Allreduce(MPI_IN_PLACE, loc_src.data(), length_src, MPI_INT, MPI_MAX, MPI_COMM_WORLD);

    for (int i = 0; i < length_src; i++)
    {
        SchurMat2DSlabX* S_src  = src != nullptr ? src->S_params[i] : nullptr;
        SchurMat2DSlabX* S_dest = nullptr;
        if (dest != nullptr)
        {
            Domain2DUniform* domain       = dest->domain;
            MPI_Comm         communicator = dest->communicator;
            switch (static_cast<LocationType>(loc_src[i]))
            {
                case LocationType::Left:
                    S_dest = new SchurMat2DSlabX_left(domain, communicator);
                    break;
                case LocationType::Right:
                    S_dest = new SchurMat2DSlabX_right(domain, communicator);
                    break;
                case LocationType::Down:
                    S_dest = new SchurMat2DSlabX_down(domain, communicator);
                    break;
                case LocationType::Up:
                    S_dest = new SchurMat2DSlabX_up(domain, communicator);
                    break;
                default:
                    std::cerr << "Unkown LocationType in GMRESSolver2DSlabX::migrate_from" << std::endl;
                    return;
                    break;
            }
        }

        migrate_from(S_src, S_dest);

        if (dest != nullptr)
            dest->S_params.push_back(S_dest);
    }
}
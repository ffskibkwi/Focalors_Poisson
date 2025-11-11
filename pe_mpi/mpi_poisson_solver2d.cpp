#include "mpi_poisson_solver2d.h"

#include <cstring>

/**
 * @brief Construct an MPI Poisson solver (直接指定网格与边界类型）。
 *
 * 与 OpenMP 版本保持数值流程一致：y-向 FFT → 全局转置 → x-向追赶法 → 反转置 → y-向逆变换。
 * 该构造函数不依赖 Variable/Domain 的边界映射，直接通过参数指定四面边界类型。
 *
 * @param in_nx    全局 x 方向尺寸。
 * @param in_ny    全局 y 方向尺寸。
 * @param in_hx    网格步长 hx。
 * @param in_hy    网格步长 hy。
 * @param in_boundary_type_left   左边界类型。
 * @param in_boundary_type_right  右边界类型。
 * @param in_boundary_type_down   下边界类型。
 * @param in_boundary_type_up     上边界类型。
 * @param in_num_proc             参与并行的进程数上限（保留参数，实际参与由外部通信器决定）。
 * @param in_start_rank           起始 rank（保留参数，不在类内进行 Comm_split）。
 * @param in_env_config           环境配置（可控制打印等）。
 * @param in_comm                 外部提供的合法 MPI 通信器（本类不再内部 Comm_split）。
 */

MPIPoissonSolver2D::MPIPoissonSolver2D(int                in_nx,
                                       int                in_ny,
                                       double             in_hx,
                                       double             in_hy,
                                       PDEBoundaryType    in_boundary_type_left,
                                       PDEBoundaryType    in_boundary_type_right,
                                       PDEBoundaryType    in_boundary_type_down,
                                       PDEBoundaryType    in_boundary_type_up,
                                       int                in_num_proc,
                                       int                in_start_rank,
                                       EnvironmentConfig* in_env_config,
                                       MPI_Comm           in_comm)
    : nx(in_nx)
    , ny(in_ny)
    , hx(in_hx)
    , hy(in_hy)
    , boundary_type_left(in_boundary_type_left)
    , boundary_type_right(in_boundary_type_right)
    , boundary_type_down(in_boundary_type_down)
    , boundary_type_up(in_boundary_type_up)
    , num_proc(in_num_proc)
    , start_rank(in_start_rank)
{
    env_config = in_env_config;
    setup_comm(in_comm);
    setup_problem_params();
    build_decomposition();
    build_local_solvers();
}

/**
 * @brief 基于 Domain/Variable 构造 MPI Poisson 求解器。
 *
 * 边界类型从 Variable 的 boundary_type_map[domain] 读取，以与当前库数据结构对齐。
 *
 * @param in_domain     计算域（提供 nx, ny, hx, hy）。
 * @param in_variable   变量（提供四面边界类型与值）。
 * @param in_num_proc   参与并行的进程数上限（保留参数，实际参与由外部通信器决定）。
 * @param in_start_rank 起始 rank（保留参数，不在类内进行 Comm_split）。
 * @param in_env_config 环境配置。
 * @param in_comm       外部提供的合法 MPI 通信器（本类不再内部 Comm_split）。
 */
MPIPoissonSolver2D::MPIPoissonSolver2D(Domain2DUniform* in_domain,
                                       Variable*        in_variable,
                                       int              in_num_proc,
                                       int              in_start_rank,
                                       EnvironmentConfig* in_env_config,
                                       MPI_Comm           in_comm)
    : domain(in_domain)
    , var(in_variable)
    , nx(in_domain->nx)
    , ny(in_domain->ny)
    , hx(in_domain->hx)
    , hy(in_domain->hy)
    , num_proc(in_num_proc)
    , start_rank(in_start_rank)
{
    env_config           = in_env_config;
    boundary_type_left   = var->boundary_type_map[domain][LocationType::Left];
    boundary_type_right  = var->boundary_type_map[domain][LocationType::Right];
    boundary_type_down   = var->boundary_type_map[domain][LocationType::Down];
    boundary_type_up     = var->boundary_type_map[domain][LocationType::Up];

    setup_comm(in_comm);
    setup_problem_params();
    build_decomposition();
    build_local_solvers();
}

/**
 * @brief 析构：释放局部求解器与子通信器。
 */
MPIPoissonSolver2D::~MPIPoissonSolver2D()
{
    release_local_solvers();
    // 外部提供的通信器不在此处释放
}

/**
 * @brief 初始化通信环境（不进行 Comm_split，直接使用外部提供的通信器）。
 *
 * 假定外部已经根据需要准备好通信子集；本类仅记录该通信器与其 rank/size。
 *
 * @param in_comm 外部传入的通信器。
 */
void MPIPoissonSolver2D::setup_comm(MPI_Comm in_comm)
{
    world_comm  = in_comm;
    active_comm = in_comm;
    if (active_comm != MPI_COMM_NULL)
    {
        is_active = true;
        MPI_Comm_rank(active_comm, &active_rank);
        MPI_Comm_size(active_comm, &active_size);
    }
    else
    {
        is_active = false;
        active_rank = -1;
        active_size = 0;
    }
}

/**
 * @brief 预留入口（当前无额外问题参数设置）。
 */
void MPIPoissonSolver2D::setup_problem_params()
{
    // nothing else for now
}

/**
 * @brief 构建一维条带分解信息。
 *
 * - i-slab：沿 x 方向（行）切分，尺寸为 (local_i_count, ny)。用于本地 y-向 FFT。
 * - j-slab：沿 y 方向（列）切分，尺寸为 (local_j_count, nx)。用于本地 x-向追赶法（在全局转置之后）。
 */
void MPIPoissonSolver2D::build_decomposition()
{
    if (!is_active)
        return;

    split_1d(nx, active_size, i_counts, i_displs);
    split_1d(ny, active_size, j_counts, j_displs);

    local_i_count = i_counts[active_rank];
    local_j_count = j_counts[active_rank];
}

/**
 * @brief 为本进程构建本地求解器（y-向 FFT 与 x-向追赶法）。
 *
 * - y 方向 FFT 的本地尺寸为 (local_i_count, ny)。
 * - x 方向追赶法在转置域进行：尺寸 (local_j_count, nx)，并基于对应全局 j 段构造 `local_x_diag`。
 */
void MPIPoissonSolver2D::build_local_solvers()
{
    if (!is_active)
        return;

    poisson_fft_y = nullptr;
    auto isDirLike = [](PDEBoundaryType t) {
        return t == PDEBoundaryType::Dirichlet || t == PDEBoundaryType::Adjacented;
    };

    // 当本地工作量为 0 时，不构建本地算子（但仍需参与 MPI 通信）
    if (local_i_count > 0 && ny > 0)
    {
        if (boundary_type_down == PDEBoundaryType::Periodic && boundary_type_up == PDEBoundaryType::Periodic)
        {
            poisson_fft_y = new PoissonFFT2D_PP();
        }
        else if (boundary_type_down == PDEBoundaryType::Neumann && boundary_type_up == PDEBoundaryType::Neumann)
        {
            poisson_fft_y = new PoissonFFT2D_NN();
        }
        else if (isDirLike(boundary_type_down) && isDirLike(boundary_type_up))
        {
            poisson_fft_y = new PoissonFFT2D_DD();
        }
        else if ((isDirLike(boundary_type_down) && boundary_type_up == PDEBoundaryType::Neumann) ||
                 (boundary_type_down == PDEBoundaryType::Neumann && isDirLike(boundary_type_up)))
        {
            poisson_fft_y = new PoissonFFT2D_DN();
        }
        if (poisson_fft_y)
            poisson_fft_y->init(local_i_count, ny);
    }

    if (local_j_count > 0)
        build_local_x_diag();

    bool is_no_Dirichlet = !(isDirLike(boundary_type_left) || isDirLike(boundary_type_right) ||
                             isDirLike(boundary_type_down) || isDirLike(boundary_type_up));
    bool has_last_vector = true;

    chasing_method_x = nullptr;
    if (local_j_count > 0 && nx > 0)
    {
        chasing_method_x = new ChasingMethod2D();
        // For transposed layout: (local_j_count, nx)
        chasing_method_x->init(local_j_count, nx, local_x_diag, is_no_Dirichlet, has_last_vector,
                               boundary_type_left, boundary_type_right);
    }

    // Allocate persistent work fields
    f_local.init(local_i_count, ny, "f_local");
    fhat_local.init(local_i_count, ny, "fhat_local");
    fhat_local_T.init(local_j_count, nx, "fhat_local_T");
    phat_local_T.init(local_j_count, nx, "phat_local_T");
    phat_local.init(local_i_count, ny, "phat_local");
    p_local.init(local_i_count, ny, "p_local");
}

/**
 * @brief 释放本地求解器与缓存。
 */
void MPIPoissonSolver2D::release_local_solvers()
{
    if (local_x_diag)
    {
        delete[] local_x_diag;
        local_x_diag = nullptr;
    }
    if (poisson_fft_y)
    {
        delete poisson_fft_y;
        poisson_fft_y = nullptr;
    }
    if (chasing_method_x)
    {
        delete chasing_method_x;
        chasing_method_x = nullptr;
    }
}

/**
 * @brief 计算 y 向特征值 lambda_y（与 OpenMP 版本一致）。
 *
 * 根据上下边界类型（Periodic/Neumann/Dirichlet/Adjacented 的组合）分别使用对应的解析公式。
 * @param lambda_y_out 输出：长度为 ny 的向量。
 */
void MPIPoissonSolver2D::cal_lambda(std::vector<double>& lambda_y_out) const
{
    lambda_y_out.resize(ny);
    for (int i = 1; i <= ny; i++)
    {
        if (boundary_type_down == PDEBoundaryType::Periodic && boundary_type_up == PDEBoundaryType::Periodic)
        {
            lambda_y_out[i - 1] = -2.0 + 2.0 * std::cos(2.0 * pi / ny * std::floor(i / 2.0));
        }
        else if (boundary_type_down == PDEBoundaryType::Neumann && boundary_type_up == PDEBoundaryType::Neumann)
        {
            lambda_y_out[i - 1] = -2.0 + 2.0 * std::cos(pi / ny * i);
        }
        else if ((boundary_type_down == PDEBoundaryType::Dirichlet || boundary_type_down == PDEBoundaryType::Adjacented) &&
                 (boundary_type_up == PDEBoundaryType::Dirichlet || boundary_type_up == PDEBoundaryType::Adjacented))
        {
            lambda_y_out[i - 1] = -2.0 + 2.0 * std::cos(pi / (ny + 1) * i);
        }
        else if (((boundary_type_down == PDEBoundaryType::Dirichlet || boundary_type_down == PDEBoundaryType::Adjacented) &&
                  boundary_type_up == PDEBoundaryType::Neumann) ||
                 (boundary_type_down == PDEBoundaryType::Neumann &&
                  (boundary_type_up == PDEBoundaryType::Dirichlet || boundary_type_up == PDEBoundaryType::Adjacented)))
        {
            lambda_y_out[i - 1] = -2.0 + 2.0 * std::cos(2.0 * pi * i / (2 * ny + 1));
        }
    }
}

/**
 * @brief 构建本地 j 段的 x 向对角项切片 local_x_diag。
 *
 * 全局公式：x_diag[j] = -2.0 + (hx^2/hy^2) * lambda_y[j]。
 * 在分布式场景中，仅为本进程拥有的 j 索引区间 [j_start, j_start + local_j_count) 生成对应切片。
 */
void MPIPoissonSolver2D::build_local_x_diag()
{
    if (!is_active)
        return;

    std::vector<double> lambda_y_vec;
    cal_lambda(lambda_y_vec);

    // Build local slice for owned j indices
    local_x_diag = new double[local_j_count];
    int j_start  = j_displs[active_rank];
    for (int j = 0; j < local_j_count; ++j)
    {
        local_x_diag[j] = -2.0 + hx * hx / hy / hy * lambda_y_vec[j_start + j];
    }
}

/**
 * @brief 在全局 root（active_rank==0）上对右端项 f 进行边界装配。
 *
 * 该步骤与 OpenMP 版本一一对应：根据 Dirichlet/Neumann 边界值修正 f 的边缘条目。
 * @param f 全局右端项（仅 root 持有完整数据）。
 */
void MPIPoissonSolver2D::boundary_assembly(field2& f)
{
    auto& var_has_map   = var->has_boundary_value_map[domain];
    auto& var_value_map = var->boundary_value_map[domain];

    if (boundary_type_left == PDEBoundaryType::Dirichlet && var_has_map[LocationType::Left])
    {
        double* boundary_value = var_value_map[LocationType::Left];
        for (int j = 0; j < ny; j++)
            f(0, j) -= boundary_value[j] / hx / hx;
    }
    if (boundary_type_right == PDEBoundaryType::Dirichlet && var_has_map[LocationType::Right])
    {
        double* boundary_value = var_value_map[LocationType::Right];
        for (int j = 0; j < ny; j++)
            f(nx - 1, j) -= boundary_value[j] / hx / hx;
    }

    if (boundary_type_down == PDEBoundaryType::Dirichlet && var_has_map[LocationType::Down])
    {
        double* boundary_value = var_value_map[LocationType::Down];
        for (int i = 0; i < nx; i++)
            f(i, 0) -= boundary_value[i] / hy / hy;
    }
    if (boundary_type_up == PDEBoundaryType::Dirichlet && var_has_map[LocationType::Up])
    {
        double* boundary_value = var_value_map[LocationType::Up];
        for (int i = 0; i < nx; i++)
            f(i, ny - 1) -= boundary_value[i] / hy / hy;
    }

    if (boundary_type_left == PDEBoundaryType::Neumann && var_has_map[LocationType::Left])
    {
        double* boundary_value = var_value_map[LocationType::Left];
        for (int j = 0; j < ny; j++)
            f(0, j) += boundary_value[j] / hx;
    }
    if (boundary_type_right == PDEBoundaryType::Neumann && var_has_map[LocationType::Right])
    {
        double* boundary_value = var_value_map[LocationType::Right];
        for (int j = 0; j < ny; j++)
            f(nx - 1, j) -= boundary_value[j] / hx;
    }
    if (boundary_type_down == PDEBoundaryType::Neumann && var_has_map[LocationType::Down])
    {
        double* boundary_value = var_value_map[LocationType::Down];
        for (int i = 0; i < nx; i++)
            f(i, 0) += boundary_value[i] / hy;
    }
    if (boundary_type_up == PDEBoundaryType::Neumann && var_has_map[LocationType::Up])
    {
        double* boundary_value = var_value_map[LocationType::Up];
        for (int i = 0; i < nx; i++)
            f(i, ny - 1) -= boundary_value[i] / hy;
    }
}

/**
 * @brief 将 root 上的全局场按 i-slab 切分并分发到各进程。
 *
 * 使用 MPI_Scatterv，单位为 double 元素数；利用 i_counts/i_displs 保证每个进程拿到连续的行块。
 * @param f_global  仅 root 持有的全局场 (nx, ny)。
 * @param f_local   输出：本地场 (local_i_count, ny)。
 */
void MPIPoissonSolver2D::scatter_global_f_to_local(const field2& f_global, field2& f_local) const
{
    // Prepare counts/displs in elements
    std::vector<int> sendcounts(active_size), displs(active_size);
    for (int r = 0; r < active_size; ++r)
    {
        sendcounts[r] = i_counts[r] * ny;
        displs[r]     = i_displs[r] * ny;
    }

    if (env_config && env_config->showCurrentStep)
    {
        if (active_rank == 0)
        {
            std::cout << "[Poisson][MPI] Scatterv: start (nx=" << nx << ", ny=" << ny << ")\n";
            for (int r = 0; r < active_size; ++r)
                std::cout << "  - rank " << r << ": sendcount=" << sendcounts[r] << " (i_count=" << i_counts[r] << ")\n";
        }
    }
    MPI_Scatterv(
        f_global.value,
        sendcounts.data(),
        displs.data(),
        MPI_DOUBLE,
        f_local.value,
        sendcounts[active_rank],
        MPI_DOUBLE,
        0,
        active_comm);
    if (env_config && env_config->showCurrentStep)
    {
        std::cout << "[Poisson][MPI] Scatterv: done at rank " << active_rank
                  << " (recv=" << sendcounts[active_rank] << ")\n";
    }
}

/**
 * @brief 将各进程的 i-slab 本地结果回收至 root 的全局场。
 *
 * 使用 MPI_Gatherv，单位为 double 元素数；与分发阶段完全对称。
 * @param f_local  本地场 (local_i_count, ny)。
 * @param f_global 仅 root 持有的全局场 (nx, ny)。
 */
void MPIPoissonSolver2D::gather_local_to_global_f(const field2& f_local, field2& f_global) const
{
    std::vector<int> recvcounts(active_size), displs(active_size);
    for (int r = 0; r < active_size; ++r)
    {
        recvcounts[r] = i_counts[r] * ny;
        displs[r]     = i_displs[r] * ny;
    }

    if (env_config && env_config->showCurrentStep)
    {
        if (active_rank == 0)
        {
            std::cout << "[Poisson][MPI] Gatherv: start\n";
            for (int r = 0; r < active_size; ++r)
                std::cout << "  - rank " << r << ": recvcount=" << recvcounts[r] << "\n";
        }
    }
    MPI_Gatherv(
        f_local.value,
        recvcounts[active_rank],
        MPI_DOUBLE,
        f_global.value,
        recvcounts.data(),
        displs.data(),
        MPI_DOUBLE,
        0,
        active_comm);
    if (env_config && env_config->showCurrentStep)
    {
        std::cout << "[Poisson][MPI] Gatherv: done at rank " << active_rank
                  << " (send=" << recvcounts[active_rank] << ")\n";
    }
}

/**
 * @brief 分布式全局转置：从 i-slab 分布 (local_i_count, ny) 转换到 j-slab 分布 (local_j_count, nx)。
 *
 * 算法说明：
 * - 按目标进程的 j 段打包每个本地 i 行的对应列段，形成按目的进程分组的 sendbuf；
 * - 通过 MPI_Alltoallv 交换数据块；
 * - 在接收端按来源进程的 i 段与本地 j 段重组为 (local_j_count, nx) 的转置布局。
 * 保证与单机 `field2::transpose` 一致的全局重新排列。
 *
 * @param in_i_slab  输入：本地 i-slab 场。
 * @param out_j_slab 输出：本地 j-slab 场（已转置）。
 */
void MPIPoissonSolver2D::distributed_transpose_i_to_j(const field2& in_i_slab, field2& out_j_slab) const
{
    // in_i_slab: (local_i_count, ny)
    // out_j_slab: (local_j_count, nx)

    // Pack send buffers per-destination rank
    std::vector<int> sendcounts(active_size), sdispls(active_size);
    std::vector<int> recvcounts(active_size), rdispls(active_size);
    for (int r = 0; r < active_size; ++r)
    {
        sendcounts[r] = local_i_count * j_counts[r];
        recvcounts[r] = i_counts[r] * local_j_count;
    }
    sdispls[0] = 0;
    rdispls[0] = 0;
    for (int r = 1; r < active_size; ++r)
    {
        sdispls[r] = sdispls[r - 1] + sendcounts[r - 1];
        rdispls[r] = rdispls[r - 1] + recvcounts[r - 1];
    }

    std::vector<double> sendbuf(sdispls.back() + sendcounts.back());
    std::vector<double> recvbuf(rdispls.back() + recvcounts.back());

    if (env_config && env_config->showCurrentStep)
    {
        std::cout << "[Poisson][MPI] Alltoallv (i->j): start at rank " << active_rank
                  << " (local_i=" << local_i_count << ", local_j=" << local_j_count << ")\n";
    }
    for (int r = 0; r < active_size; ++r)
    {
        int j0 = j_displs[r];
        int jc = j_counts[r];
        double* dst = sendbuf.data() + sdispls[r];
        for (int i = 0; i < local_i_count; ++i)
        {
            // 每一行 i，复制该行中属于目标进程 r 的 j 段 [j0, j0 + jc)
            std::memcpy(dst + i * jc, in_i_slab.get_ptr(i, j0), sizeof(double) * jc);
        }
    }

    MPI_Alltoallv(sendbuf.data(), sendcounts.data(), sdispls.data(), MPI_DOUBLE,
                  recvbuf.data(), recvcounts.data(), rdispls.data(), MPI_DOUBLE,
                  active_comm);
    if (env_config && env_config->showCurrentStep)
    {
        std::cout << "[Poisson][MPI] Alltoallv (i->j): done at rank " << active_rank << "\n";
    }

    // Unpack into out_j_slab (local_j_count, nx)
    out_j_slab.set_size(local_j_count, nx);
    for (int s = 0; s < active_size; ++s)
    {
        int ic = i_counts[s];
        int i0 = i_displs[s];
        const double* src = recvbuf.data() + rdispls[s];
        for (int il = 0; il < ic; ++il)
        {
            int i_global = i0 + il;
            for (int jl = 0; jl < local_j_count; ++jl)
            {
                double v = src[il * local_j_count + jl];
                out_j_slab(jl, i_global) = v;
            }
        }
    }
}

/**
 * @brief 分布式全局反转置：从 j-slab 分布 (local_j_count, nx) 回到 i-slab 分布 (local_i_count, ny)。
 *
 * 与 distributed_transpose_i_to_j 对称：
 * - 发送端按每个目标进程的 i 段，从 (local_j_count, nx) 中抽取整列块；
 * - MPI_Alltoallv 交换；
 * - 接收端按 (local_i_count, ny) 重组。
 *
 * @param in_j_slab  输入：本地 j-slab 场。
 * @param out_i_slab 输出：本地 i-slab 场。
 */
void MPIPoissonSolver2D::distributed_transpose_j_to_i(const field2& in_j_slab, field2& out_i_slab) const
{
    // in_j_slab: (local_j_count, nx)
    // out_i_slab: (local_i_count, ny)

    std::vector<int> sendcounts(active_size), sdispls(active_size);
    std::vector<int> recvcounts(active_size), rdispls(active_size);
    for (int r = 0; r < active_size; ++r)
    {
        sendcounts[r] = i_counts[r] * local_j_count;
        recvcounts[r] = local_i_count * j_counts[r];
    }
    sdispls[0] = 0;
    rdispls[0] = 0;
    for (int r = 1; r < active_size; ++r)
    {
        sdispls[r] = sdispls[r - 1] + sendcounts[r - 1];
        rdispls[r] = rdispls[r - 1] + recvcounts[r - 1];
    }

    std::vector<double> sendbuf(sdispls.back() + sendcounts.back());
    std::vector<double> recvbuf(rdispls.back() + recvcounts.back());

    if (env_config && env_config->showCurrentStep)
    {
        std::cout << "[Poisson][MPI] Alltoallv (j->i): start at rank " << active_rank
                  << " (local_i=" << local_i_count << ", local_j=" << local_j_count << ")\n";
    }
    for (int r = 0; r < active_size; ++r)
    {
        int i0 = i_displs[r];
        int ic = i_counts[r];
        double* dst = sendbuf.data() + sdispls[r];
        for (int il = 0; il < ic; ++il)
        {
            int i_global = i0 + il;
            for (int jl = 0; jl < local_j_count; ++jl)
            {
                // 取列 i_global 的本地 j 段，线性化为 (ic x local_j_count) 的块
                dst[il * local_j_count + jl] = in_j_slab(jl, i_global);
            }
        }
    }

    MPI_Alltoallv(sendbuf.data(), sendcounts.data(), sdispls.data(), MPI_DOUBLE,
                  recvbuf.data(), recvcounts.data(), rdispls.data(), MPI_DOUBLE,
                  active_comm);
    if (env_config && env_config->showCurrentStep)
    {
        std::cout << "[Poisson][MPI] Alltoallv (j->i): done at rank " << active_rank << "\n";
    }

    out_i_slab.set_size(local_i_count, ny);
    for (int s = 0; s < active_size; ++s)
    {
        int jc = j_counts[s];
        int j0 = j_displs[s];
        const double* src = recvbuf.data() + rdispls[s];
        for (int il = 0; il < local_i_count; ++il)
        {
            for (int jl = 0; jl < jc; ++jl)
            {
                // 将 (ic x jc) 的块按 j 全局偏移 j0 放入 out_i_slab 的列段
                out_i_slab(il, j0 + jl) = src[il * jc + jl];
            }
        }
    }
}

/**
 * @brief 求解主流程（仅活动通信器进程参与）。
 *
 * 流程保持与单机 OpenMP 版本一致：
 * 1) root 做边界装配；2) 按 i-slab 分发；3) 本地 y-FFT；4) 全局转置 i→j；
 * 5) 本地 x-追赶；6) 全局转置 j→i；7) 本地 y-逆变换；8) 回收至 root。
 * @param f 仅 root 持有的全局场（输入右端项，输出解）。
 */
void MPIPoissonSolver2D::solve(field2& f)
{
    // Only active ranks participate
    if (!is_active)
        return;

    if (env_config && env_config->showCurrentStep && active_rank == 0)
    {
        std::cout << "[Poisson][MPI] solve: start" << std::endl;
        double s_in = f.sum();
        std::cout << "[Poisson][MPI] root input f.sum(before asm)=" << s_in << std::endl;
    }

    // 1) Boundary assembly on root, since full field is only on rank 0
    if (active_rank == 0 && var && domain)
    {
        boundary_assembly(f);
        if (env_config && env_config->showCurrentStep)
        {
            double s_after = f.sum();
            std::cout << "[Poisson][MPI] root input f.sum(after asm)=" << s_after << std::endl;
        }
    }

    if (env_config && env_config->showCurrentStep)
    {
        std::cout << "[Poisson][MPI] rank " << active_rank
                  << " local_i=" << local_i_count << " local_j=" << local_j_count
                  << " nx=" << nx << " ny=" << ny << std::endl;
    }

    // 2) Scatter to i-slab local fields
    scatter_global_f_to_local(f, f_local);

    // 3) Y-direction transform locally: f_local -> fhat_local
    if (poisson_fft_y && local_i_count > 0)
        poisson_fft_y->transform(f_local, fhat_local);

    // 4) Global transpose: i-slab -> j-slab
    distributed_transpose_i_to_j(fhat_local, fhat_local_T);

    // 5) X-direction chasing on transposed slabs (local_j_count, nx)
    if (chasing_method_x && local_j_count > 0)
        chasing_method_x->chasing(fhat_local_T, phat_local_T);

    // 6) Global transpose back: j-slab -> i-slab
    distributed_transpose_j_to_i(phat_local_T, phat_local);

    // 7) Inverse Y-transform locally
    if (poisson_fft_y && local_i_count > 0)
        poisson_fft_y->transform_transpose(phat_local, p_local);

    // 8) Gather back to root
    if (env_config && env_config->showCurrentStep && active_rank == 0)
    {
        double s_local = p_local.sum();
        std::cout << "[Poisson][MPI] rank0 pre-gather p_local.sum=" << s_local
                  << " (nx=" << nx << ", ny=" << ny << ", li=" << local_i_count << ", lj=" << local_j_count << ")\n";
    }
    gather_local_to_global_f(p_local, f);

    if (env_config && env_config->showCurrentStep && active_rank == 0)
    {
        double s_f = f.sum();
        std::cout << "[Poisson][MPI] solve: done, f.sum=" << s_f << std::endl;
    }
}

void MPIPoissonSolver2D::solve_collective_root_owned(field2& f)
{
    if (!is_active)
        return;
    if (active_rank == 0)
    {
        solve(f);
    }
    else
    {
        // 非 root：提供占位缓冲参与通信
        field2 dummy;
        dummy.init(nx, ny);
        solve(dummy);
    }
}

/**
 * @brief 将长度为 n 的一维区间均匀切分给 p 个进程，生成 counts 与 displs。
 *
 * @param n       全局长度。
 * @param p       进程数。
 * @param counts  输出：每个进程的长度。
 * @param displs  输出：每个进程的起始偏移。
 */
void MPIPoissonSolver2D::split_1d(int n, int p, std::vector<int>& counts, std::vector<int>& displs)
{
    counts.assign(p, 0);
    displs.assign(p, 0);
    int base = n / p;
    int rem  = n % p;
    int off  = 0;
    for (int r = 0; r < p; ++r)
    {
        counts[r] = base + (r < rem ? 1 : 0);
        displs[r] = off;
        off += counts[r];
    }
}



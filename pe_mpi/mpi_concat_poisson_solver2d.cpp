#include "mpi_concat_poisson_solver2d.h"

#include <algorithm>
#include <stdexcept>
#include <iostream>
#include <sstream>

MPIConcatPoissonSolver2D::MPIConcatPoissonSolver2D(Variable* in_variable, MPI_Comm in_comm, EnvironmentConfig* in_env_config)
    : variable(in_variable)
    , env_config(in_env_config)
    , global_comm(in_comm)
{
    if (variable == nullptr)
        throw std::runtime_error("MPIConcatPoissonSolver2D: variable is null");
    if (variable->geometry == nullptr)
        throw std::runtime_error("MPIConcatPoissonSolver2D: variable->geometry is null");

    geometry  = variable->geometry;
    // 确保 geometry 已检查与准备（避免 Geometry2D is not checked）
    if (!geometry->is_checked)
        geometry->check();
    if (!geometry->is_prepared)
        geometry->solve_prepare();

    tree_root  = geometry->tree_root;
    tree_map   = geometry->tree_map;
    parent_map = geometry->parent_map;
    field_map  = variable->field_map;

    MPI_Comm_rank(global_comm, &world_rank);
    MPI_Comm_size(global_comm, &world_size);

    build_levels();
    build_domain_communicators();
    construct_solver_map();

    // 构建临时场（与串行一致，非根域使用临时 b）
    for (auto &[domain, fptr] : field_map)
    {
        if (domain != tree_root)
            temp_fields[domain] = new field2(fptr->get_nx(), fptr->get_ny(), fptr->get_name() + "_temp");
    }
}

MPIConcatPoissonSolver2D::~MPIConcatPoissonSolver2D()
{
    for (auto &[domain, temp_field] : temp_fields)
        delete temp_field;
    for (auto &domain : field_map)
    {
        if (solver_map[domain.first])
            delete solver_map[domain.first];
    }
    for (auto &[domain, comm] : domain_comm_map)
    {
        if (comm != MPI_COMM_NULL && comm != global_comm)
            MPI_Comm_free(&comm);
    }
}

/**
 * @brief 层级构建：BFS，root 层为 0。
 */
void MPIConcatPoissonSolver2D::build_levels()
{
    levels = TreeUtils::buildLevelsFromTree(tree_root, tree_map);
}

/**
 * @brief 将世界进程按层进行平均划分，并为每个域创建子通信器。
 *
 * 层 0（root 域）使用全部进程；层 l>0 上有 k 个域时，将 world_size 平均分为 k 份（向前取整），
 * 按连续 rank 区间为每个域创建一个子通信器。
 */
void MPIConcatPoissonSolver2D::build_domain_communicators()
{
    // 关键修正：通信器按树形层级“父域范围内嵌套”划分，确保子域通信器完全包含于父域通信器
    // 层 0：root 使用全体进程
    domain_comm_map[tree_root] = global_comm;
    domain_world_range[tree_root] = std::make_pair(0, world_size);

    // 层 >0：对每个父域的已分配进程区间进行二次划分给其子域
    for (size_t lev = 1; lev < levels.size(); ++lev)
    {
        // 收集本层各父域的子域列表
        std::unordered_map<Domain2DUniform*, std::vector<Domain2DUniform*>> parent_to_children;
        for (auto* d : levels[lev])
        {
            auto itp = parent_map.find(d);
            if (itp == parent_map.end()) continue;
            Domain2DUniform* parent = itp->second.second;
            parent_to_children[parent].push_back(d);
        }

        // 针对每个父域，在其已分配范围内均分给子域
        for (auto &kv : parent_to_children)
        {
            Domain2DUniform* parent = kv.first;
            auto &children = kv.second;
            if (domain_world_range.find(parent) == domain_world_range.end())
                continue; // 安全保护

            int p_start = domain_world_range[parent].first;
            int p_cnt   = domain_world_range[parent].second;

            std::vector<int> counts, displs;
            even_split_ranges(p_cnt, static_cast<int>(children.size()), counts, displs);

            for (size_t idx = 0; idx < children.size(); ++idx)
            {
                Domain2DUniform* child = children[idx];
                int c_start = p_start + displs[idx];
                int c_cnt   = counts[idx];

                // 使用“子域起始 world rank”作为 color，保证全局唯一
                int color = (world_rank >= c_start && world_rank < c_start + c_cnt) ? c_start : MPI_UNDEFINED;
                MPI_Comm sub = MPI_COMM_NULL;
                MPI_Comm_split(global_comm, color, world_rank, &sub);
                domain_comm_map[child] = sub;
                domain_world_range[child] = std::make_pair(c_start, c_cnt);
            }
        }
    }
    // 打印规划
    print_parallel_plan();
}

void MPIConcatPoissonSolver2D::even_split_ranges(int nProc, int nParts, std::vector<int>& counts, std::vector<int>& displs)
{
    counts.assign(nParts, 0);
    displs.assign(nParts, 0);
    int base = nProc / nParts;
    int rem  = nProc % nParts;
    int off  = 0;
    for (int i = 0; i < nParts; ++i)
    {
        counts[i] = base + (i < rem ? 1 : 0);
        displs[i] = off;
        off += counts[i];
    }
}

/**
 * @brief 构造 MPI 版 solver_map：叶子→MPIPoisson；非叶→MPIGMRES（内部也用 MPI Poisson）。
 */
void MPIConcatPoissonSolver2D::construct_solver_map()
{
    // 关键修正：从深层到浅层构建，使子域 solver 先于父域可用（与串行一致）
    for (int lev = static_cast<int>(levels.size()) - 1; lev >= 0; --lev)
    {
        for (auto *domain : levels[static_cast<size_t>(lev)])
        {
            MPI_Comm dcomm = domain_comm_map[domain];
            if (dcomm == MPI_COMM_NULL)
            {
                // 本 rank 不参与该域的求解器通信器，避免在 NULL 通信器上调用 MPI 接口
                solver_map[domain] = nullptr;
                continue;
            }
            if (tree_map[domain].size() > 0) // 非叶
            {
                // GMRES（MPI）
                // 数值保持与串行一致：m=20, tol=1e-8, maxIter=100
                solver_map[domain] = new MPIGMRESSolver2D(domain, variable, /*m*/20, /*tol*/1e-8, /*maxIter*/100, dcomm, env_config);
            }
            else
            {
                // 叶子直接 Poisson（MPI）
                int dsize = 0;
                MPI_Comm_size(dcomm, &dsize);
                auto* pe = new MPIPoissonSolver2D(domain, variable, dsize, /*start_rank*/0, env_config, dcomm);
                solver_map[domain] = pe;
            }
            if (env_config && env_config->showCurrentStep)
            {
                int rr = -1; MPI_Comm_rank(dcomm, &rr);
                if (rr == 0)
                    std::cout << "[MPIConcat] Built solver for domain " << domain->name << " in its subcomm" << std::endl;
            }

            // 在确保子域 solver 已经存在后，再构建当前非叶域的 Schur
            if (tree_map[domain].size() > 0)
            {
                auto* gmres = static_cast<MPIGMRESSolver2D*>(solver_map[domain]);
                if (gmres)
                    gmres->SchurMat2D_construct(tree_map[domain], solver_map);
            }
        }
    }
}

/**
 * @brief 获取子通信器的 local rank 0 对应的 world rank。
 */
int MPIConcatPoissonSolver2D::comm_rank_of_local_root(MPI_Comm subcomm) const
{
    if (subcomm == MPI_COMM_NULL) return -1;
    MPI_Group world_group, sub_group;
    MPI_Comm_group(global_comm, &world_group);
    MPI_Comm_group(subcomm, &sub_group);

    int root_in_sub = 0;
    int world_root  = -1;
    MPI_Group_translate_ranks(sub_group, 1, &root_in_sub, world_group, &world_root);

    MPI_Group_free(&sub_group);
    MPI_Group_free(&world_group);
    return world_root;
}

/**
 * @brief 将全局 root 上的场推送到子通信器的 root（若两者不同）。
 */
void MPIConcatPoissonSolver2D::push_field_to_comm_root(Domain2DUniform* domain, MPI_Comm subcomm, field2& f_global) const
{
    int sub_rank = -1, sub_size = 0;
    if (subcomm != MPI_COMM_NULL)
    {
        MPI_Comm_rank(subcomm, &sub_rank);
        MPI_Comm_size(subcomm, &sub_size);
    }

    // 计算子通信器 root 的 world rank：优先从规划中读取
    int sub_root_world = -1;
    auto it = domain_world_range.find(domain);
    if (it != domain_world_range.end())
        sub_root_world = it->second.first;
    else if (subcomm != MPI_COMM_NULL)
        sub_root_world = comm_rank_of_local_root(subcomm);

    if (world_rank == 0 && sub_root_world == 0)
        return; // 同一进程，无需传输

    if (world_rank == 0 && sub_root_world >= 0)
    {
        MPI_Send(f_global.value, f_global.get_size_n(), MPI_DOUBLE, sub_root_world, 101, global_comm);
    }
    if (subcomm != MPI_COMM_NULL && world_rank == sub_root_world && sub_rank == 0)
    {
        MPI_Recv(f_global.value, f_global.get_size_n(), MPI_DOUBLE, 0, 101, global_comm, MPI_STATUS_IGNORE);
    }
}

/**
 * @brief 规约到 world_rank==0 打印并行规划（每层、每域的进程区间）。
 */
void MPIConcatPoissonSolver2D::print_parallel_plan() const
{
    // 仅 root 收集与打印
    if (world_rank != 0) return;
    std::ostringstream oss;
    oss << "[MPIConcat] Parallel plan (world_size=" << world_size << ")\n";
    for (size_t lev = 0; lev < levels.size(); ++lev)
    {
        oss << "  Level " << lev << ":\n";
        for (auto* d : levels[lev])
        {
            auto it = domain_world_range.find(d);
            if (it != domain_world_range.end())
            {
                int start = it->second.first;
                int cnt   = it->second.second;
                oss << "    - Domain " << d->name << ": ranks [" << start << ", " << (start+cnt-1) << "] (size=" << cnt << ")\n";
            }
            else
            {
                oss << "    - Domain " << d->name << ": (no comm)\n";
            }
        }
    }
    std::cout << oss.str() << std::flush;
}

/**
 * @brief 从子通信器 root 拉回结果到全局 root（若两者不同）。
 */
void MPIConcatPoissonSolver2D::pull_field_from_comm_root(Domain2DUniform* domain, MPI_Comm subcomm, field2& f_global) const
{
    int sub_rank = -1, sub_size = 0;
    if (subcomm != MPI_COMM_NULL)
    {
        MPI_Comm_rank(subcomm, &sub_rank);
        MPI_Comm_size(subcomm, &sub_size);
    }

    int sub_root_world = -1;
    auto it = domain_world_range.find(domain);
    if (it != domain_world_range.end())
        sub_root_world = it->second.first;
    else if (subcomm != MPI_COMM_NULL)
        sub_root_world = comm_rank_of_local_root(subcomm);

    if (world_rank == 0 && sub_root_world == 0)
        return; // 同一进程，无需传输

    if (subcomm != MPI_COMM_NULL && world_rank == sub_root_world && sub_rank == 0)
    {
        MPI_Send(f_global.value, f_global.get_size_n(), MPI_DOUBLE, 0, 102, global_comm);
    }
    if (world_rank == 0 && sub_root_world >= 0)
    {
        MPI_Recv(f_global.value, f_global.get_size_n(), MPI_DOUBLE, sub_root_world, 102, global_comm, MPI_STATUS_IGNORE);
    }
}

/**
 * @brief 右端项组装阶段：与串行一致，使用 temp_fields，并对每个域调用对应的 solver。
 *
 * 按层并行：每一层内多个域各自在其子通信器上并行调用求解。数据搬移：全局 root 与域子通信器 root 之间双向传输。
 */
void MPIConcatPoissonSolver2D::solve_right_hand_construction()
{
    // 正确顺序：自底向上（保证先解出子域，再用其影响组装父域右端）
    if (levels.empty()) return;
    for (int lev = static_cast<int>(levels.size()) - 1; lev >= 1; --lev)
    {
        auto& domains = levels[static_cast<size_t>(lev)];
        if (env_config && env_config->showCurrentStep && world_rank == 0)
        {
            std::ostringstream oss;
            oss << "[MPIConcat] RHS level " << lev << ": " << domains.size() << " domain(s)\n";
            for (auto* d : domains)
            {
                auto it = domain_world_range.find(d);
                if (it != domain_world_range.end())
                    oss << "  - " << d->name << " ranks [" << it->second.first << ", " << (it->second.first + it->second.second - 1) << "]\n";
                else
                    oss << "  - " << d->name << " (no comm)\n";
            }
            std::cout << oss.str() << std::flush;
        }
        // 预先在全局 root 上构造 temp_fields[domain]
        if (world_rank == 0)
        {
            for (auto* domain : domains)
            {
                (*temp_fields[domain]) = (*field_map[domain]);
                for (auto &[location, child_domain] : tree_map[domain])
                {
                    temp_fields[domain]->bond_add(location, -1., *temp_fields[child_domain]);
                    field_map[domain]->bond_add(location, -1., *temp_fields[child_domain]);
                }
                if (env_config && env_config->showCurrentStep)
                {
                    double s_pre = temp_fields[domain]->sum();
                    std::cout << "[MPIConcat] Domain " << domain->name << " temp sum before solve=" << s_pre << std::endl;
                }
            }
        }

        // 分别在各自子通信器上进行 Poisson/GMRES 求解
        for (auto* domain : domains)
        {
            MPI_Comm dcomm = domain_comm_map[domain];
            // 将右端项在 world_root 与子通信器 root 之间传输（所有进程均进入以匹配双向P2P）
            push_field_to_comm_root(domain, dcomm, *temp_fields[domain]);
            // 在子通信器中调用求解
            if (dcomm != MPI_COMM_NULL)
            {
                int dr = -1; MPI_Comm_rank(dcomm, &dr);
                // 所有子通信器内进程均需进入 solve，以参与其内部的 MPI 通信
                if (env_config && env_config->showCurrentStep && dr == 0)
                    std::cout << "[MPIConcat] Domain " << domain->name << " solve: start" << std::endl;
                field2 dummy;
                if (dr != 0)
                    dummy.init(domain->nx, domain->ny); // 非 root 提供占位缓冲
                field2& arg = (dr == 0) ? *temp_fields[domain] : dummy;
                if (solver_map[domain])
                    solver_map[domain]->solve(arg);
                if (env_config && env_config->showCurrentStep && dr == 0)
                    std::cout << "[MPIConcat] Domain " << domain->name << " solve: done" << std::endl;
            }
            // 将解从子通信器 root 拉回全局 root（所有进程均进入以匹配双向P2P）
            pull_field_from_comm_root(domain, dcomm, *temp_fields[domain]);
            if (env_config && env_config->showCurrentStep && world_rank == 0)
            {
                double s_post = temp_fields[domain]->sum();
                std::cout << "[MPIConcat] Domain " << domain->name << " temp sum after solve=" << s_post << std::endl;
            }
        }
        if (env_config && env_config->showCurrentStep && world_rank == 0)
            std::cout << "[MPIConcat] RHS level " << lev << ": done" << std::endl;
    }
}

/**
 * @brief 根域方程：将子域影响加到 root 域右端上，并在 root 子通信器上求解。
 */
void MPIConcatPoissonSolver2D::solve_root_equation()
{
    if (env_config && env_config->showCurrentStep && world_rank == 0)
        std::cout << "[MPIConcat] Root equation: start" << std::endl;
    if (world_rank == 0)
    {
        for (auto &[location, child_domain] : tree_map[tree_root])
            field_map[tree_root]->bond_add(location, -1., *temp_fields[child_domain]);
        if (env_config && env_config->showCurrentStep)
        {
            double s_root_pre = field_map[tree_root]->sum();
            std::cout << "[MPIConcat] Root domain " << tree_root->name << " field sum before solve=" << s_root_pre << std::endl;
        }
    }
    // 推送到根通信器（根层使用所有进程进入，以匹配点对点）
    push_field_to_comm_root(tree_root, domain_comm_map[tree_root], *field_map[tree_root]);
    if (domain_comm_map[tree_root] != MPI_COMM_NULL)
    {
        int rr = -1; MPI_Comm_rank(domain_comm_map[tree_root], &rr);
        if (env_config && env_config->showCurrentStep && rr == 0)
            std::cout << "[MPIConcat] Root domain " << tree_root->name << " solve: start" << std::endl;
        field2 dummy;
        if (rr != 0)
            dummy.init(tree_root->nx, tree_root->ny);
        field2& arg = (rr == 0) ? *field_map[tree_root] : dummy;
        if (solver_map[tree_root])
            solver_map[tree_root]->solve(arg);
        if (env_config && env_config->showCurrentStep && rr == 0)
            std::cout << "[MPIConcat] Root domain " << tree_root->name << " solve: done" << std::endl;
    }
    // 拉回到 root（所有进程进入，以匹配点对点）
    pull_field_from_comm_root(tree_root, domain_comm_map[tree_root], *field_map[tree_root]);
    if (env_config && env_config->showCurrentStep && world_rank == 0)
    {
        double s_root_post = field_map[tree_root]->sum();
        std::cout << "[MPIConcat] Root domain " << tree_root->name << " field sum after solve=" << s_root_post << std::endl;
    }
    if (env_config && env_config->showCurrentStep && world_rank == 0)
        std::cout << "[MPIConcat] Root equation: done" << std::endl;
}

/**
 * @brief 支域回代阶段：自底向上回代求解每个支域。
 */
void MPIConcatPoissonSolver2D::solve_branch_equations()
{
    if (env_config && env_config->showCurrentStep && world_rank == 0)
        std::cout << "[MPIConcat] Branch equations: start" << std::endl;
    // 串行逻辑：由浅到深（父域先于子域），保证回代时父域解已就绪
    for (size_t lev = 1; lev < levels.size(); ++lev)
    {
        if (env_config && env_config->showCurrentStep && world_rank == 0)
        {
            std::ostringstream oss;
            oss << "[MPIConcat] Branch level " << lev << ": " << levels[lev].size() << " domain(s)\n";
            for (auto* d : levels[lev])
            {
                auto it = domain_world_range.find(d);
                if (it != domain_world_range.end())
                    oss << "  - " << d->name << " ranks [" << it->second.first << ", " << (it->second.first + it->second.second - 1) << "]\n";
                else
                    oss << "  - " << d->name << " (no comm)\n";
            }
            std::cout << oss.str() << std::flush;
        }
        // 同层并行
        for (auto* d : levels[lev])
        {
            if (world_rank == 0)
            {
                field_map[d]->bond_add(parent_map[d].first, -1., *field_map[parent_map[d].second]);
                if (env_config && env_config->showCurrentStep)
                {
                    double s_pre = field_map[d]->sum();
                    std::cout << "[MPIConcat] Branch domain " << d->name << " field sum before solve=" << s_pre << std::endl;
                }
            }
            push_field_to_comm_root(d, domain_comm_map[d], *field_map[d]);
            if (domain_comm_map[d] != MPI_COMM_NULL)
            {
                int dr = -1; MPI_Comm_rank(domain_comm_map[d], &dr);
                if (env_config && env_config->showCurrentStep && dr == 0)
                    std::cout << "[MPIConcat] Branch domain " << d->name << " solve: start" << std::endl;
                field2 dummy;
                if (dr != 0)
                    dummy.init(d->nx, d->ny);
                field2& arg = (dr == 0) ? *field_map[d] : dummy;
                if (solver_map[d])
                    solver_map[d]->solve(arg);
                if (env_config && env_config->showCurrentStep && dr == 0)
                    std::cout << "[MPIConcat] Branch domain " << d->name << " solve: done" << std::endl;
            }
            pull_field_from_comm_root(d, domain_comm_map[d], *field_map[d]);
            if (env_config && env_config->showCurrentStep && world_rank == 0)
            {
                double s_post = field_map[d]->sum();
                std::cout << "[MPIConcat] Branch domain " << d->name << " field sum after solve=" << s_post << std::endl;
            }
        }
        if (env_config && env_config->showCurrentStep && world_rank == 0)
            std::cout << "[MPIConcat] Branch level " << lev << ": done" << std::endl;
    }
    if (env_config && env_config->showCurrentStep && world_rank == 0)
        std::cout << "[MPIConcat] Branch equations: done" << std::endl;
}

void MPIConcatPoissonSolver2D::solve()
{
    // 与串行 Concat 流程对应的三阶段
    if (env_config && env_config->showCurrentStep && world_rank == 0)
        std::cout << "[MPIConcat] solve: start" << std::endl;
    solve_right_hand_construction();
    solve_root_equation();
    solve_branch_equations();
    if (env_config && env_config->showCurrentStep && world_rank == 0)
        std::cout << "[MPIConcat] solve: done" << std::endl;
}



#include "mpi_concat_poisson_solver2d.h"

#include <algorithm>
#include <stdexcept>
#include <iostream>

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
    tree_root = geometry->tree_root;
    if (tree_root == nullptr || geometry->tree_map.empty())
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
    std::unordered_map<Domain2DUniform*, int> level_of;
    std::queue<Domain2DUniform*> q;
    q.push(tree_root);
    level_of[tree_root] = 0;
    int maxL = 0;

    while (!q.empty())
    {
        auto cur = q.front(); q.pop();
        int lc = level_of[cur];
        maxL = std::max(maxL, lc);
        if (tree_map.count(cur))
        {
            for (auto &kv : tree_map[cur])
            {
                auto child = kv.second;
                level_of[child] = lc + 1;
                q.push(child);
            }
        }
    }
    levels.assign(maxL + 1, {});
    for (auto &kv : level_of)
        levels[kv.second].push_back(kv.first);
}

/**
 * @brief 将世界进程按层进行平均划分，并为每个域创建子通信器。
 *
 * 层 0（root 域）使用全部进程；层 l>0 上有 k 个域时，将 world_size 平均分为 k 份（向前取整），
 * 按连续 rank 区间为每个域创建一个子通信器。
 */
void MPIConcatPoissonSolver2D::build_domain_communicators()
{
    // 层 0：root 使用全体进程
    domain_comm_map[tree_root] = global_comm;

    // 层 >0：平均划分
    for (size_t lev = 1; lev < levels.size(); ++lev)
    {
        auto &domains = levels[lev];
        int k = static_cast<int>(domains.size());
        if (k <= 0) continue;

        std::vector<int> counts, displs;
        even_split_ranges(world_size, k, counts, displs);

        for (int i = 0; i < k; ++i)
        {
            int start = displs[i];
            int cnt   = counts[i];
            int color = (world_rank >= start && world_rank < start + cnt) ? i : MPI_UNDEFINED;
            MPI_Comm sub = MPI_COMM_NULL;
            MPI_Comm_split(global_comm, color, world_rank, &sub);
            domain_comm_map[domains[i]] = sub;
        }
    }
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
    for (auto &level_vec : levels)
    {
        for (auto *domain : level_vec)
        {
            MPI_Comm dcomm = domain_comm_map[domain];
            if (tree_map[domain].size() > 0) // 非叶
            {
                // GMRES（MPI）
                auto* gmres = new MPIGMRESSolver2D(domain, variable, /*m*/30, /*tol*/1e-8, /*maxIter*/50, dcomm, env_config);
                solver_map[domain] = gmres;
            }
            else
            {
                // 叶子直接 Poisson（MPI）
                int dsize = 0;
                if (dcomm != MPI_COMM_NULL) MPI_Comm_size(dcomm, &dsize);
                auto* pe = new MPIPoissonSolver2D(domain, variable, dsize, /*start_rank*/0, env_config, dcomm);
                solver_map[domain] = pe;
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
void MPIConcatPoissonSolver2D::push_field_to_comm_root(MPI_Comm subcomm, field2& f_global) const
{
    if (subcomm == MPI_COMM_NULL) return;
    int sub_rank = -1, sub_size = 0;
    MPI_Comm_rank(subcomm, &sub_rank);
    MPI_Comm_size(subcomm, &sub_size);
    if (sub_size <= 0) return;

    int sub_root_world = comm_rank_of_local_root(subcomm);
    if (world_rank == 0 && sub_root_world == 0)
        return; // 同一进程，无需传输

    if (world_rank == 0 && sub_root_world >= 0)
    {
        MPI_Send(f_global.value, f_global.get_size_n(), MPI_DOUBLE, sub_root_world, 101, global_comm);
    }
    if (world_rank == sub_root_world && sub_rank == 0)
    {
        MPI_Recv(f_global.value, f_global.get_size_n(), MPI_DOUBLE, 0, 101, global_comm, MPI_STATUS_IGNORE);
    }
}

/**
 * @brief 从子通信器 root 拉回结果到全局 root（若两者不同）。
 */
void MPIConcatPoissonSolver2D::pull_field_from_comm_root(MPI_Comm subcomm, field2& f_global) const
{
    if (subcomm == MPI_COMM_NULL) return;
    int sub_rank = -1, sub_size = 0;
    MPI_Comm_rank(subcomm, &sub_rank);
    MPI_Comm_size(subcomm, &sub_size);
    if (sub_size <= 0) return;

    int sub_root_world = comm_rank_of_local_root(subcomm);
    if (world_rank == 0 && sub_root_world == 0)
        return; // 同一进程，无需传输

    if (world_rank == sub_root_world && sub_rank == 0)
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
    for (size_t lev = 1; lev < levels.size(); ++lev) // 非 root 层
    {
        auto& domains = levels[lev];
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
            }
        }

        // 分别在各自子通信器上进行 Poisson/GMRES 求解
        for (auto* domain : domains)
        {
            MPI_Comm dcomm = domain_comm_map[domain];
            // 将右端项发送到子通信器 root
            if (world_rank == 0)
                push_field_to_comm_root(dcomm, *temp_fields[domain]);
            // 在子通信器中调用求解
            if (dcomm != MPI_COMM_NULL)
            {
                int dr = -1; MPI_Comm_rank(dcomm, &dr);
                if (dr == 0)
                {
                    solver_map[domain]->solve(*temp_fields[domain]);
                }
            }
            // 将解拉回全局 root
            if (world_rank == 0)
                pull_field_from_comm_root(dcomm, *temp_fields[domain]);
        }
    }
}

/**
 * @brief 根域方程：将子域影响加到 root 域右端上，并在 root 子通信器上求解。
 */
void MPIConcatPoissonSolver2D::solve_root_equation()
{
    if (world_rank == 0)
    {
        for (auto &[location, child_domain] : tree_map[tree_root])
            field_map[tree_root]->bond_add(location, -1., *temp_fields[child_domain]);
        // 推送到根通信器（根层使用所有进程）
        push_field_to_comm_root(domain_comm_map[tree_root], *field_map[tree_root]);
    }
    if (domain_comm_map[tree_root] != MPI_COMM_NULL)
    {
        int rr = -1; MPI_Comm_rank(domain_comm_map[tree_root], &rr);
        if (rr == 0)
            solver_map[tree_root]->solve(*field_map[tree_root]);
    }
    if (world_rank == 0)
    {
        pull_field_from_comm_root(domain_comm_map[tree_root], *field_map[tree_root]);
    }
}

/**
 * @brief 支域回代阶段：自底向上回代求解每个支域。
 */
void MPIConcatPoissonSolver2D::solve_branch_equations()
{
    for (int lev = static_cast<int>(levels.size()) - 1; lev >= 1; --lev)
    {
        // 同层并行
        for (auto* d : levels[lev])
        {
            if (world_rank == 0)
            {
                field_map[d]->bond_add(parent_map[d].first, -1., *field_map[parent_map[d].second]);
                push_field_to_comm_root(domain_comm_map[d], *field_map[d]);
            }
            if (domain_comm_map[d] != MPI_COMM_NULL)
            {
                int dr = -1; MPI_Comm_rank(domain_comm_map[d], &dr);
                if (dr == 0)
                    solver_map[d]->solve(*field_map[d]);
            }
            if (world_rank == 0)
                pull_field_from_comm_root(domain_comm_map[d], *field_map[d]);
        }
    }
}

void MPIConcatPoissonSolver2D::solve()
{
    // 与串行 Concat 流程对应的三阶段
    solve_right_hand_construction();
    solve_root_equation();
    solve_branch_equations();
}



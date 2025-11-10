#pragma once

#include "base/pch.h"

#include "pe/concat/domain_solver.h"
#include "base/domain/geometry2d.h"
#include "base/domain/domain2d.h"
#include "base/domain/variable.h"
#include "io/config.h"

#include "mpi_poisson_solver2d.h"
#include "mpi_gmres_solver2d.h"

#include <mpi.h>
#include <queue>
#include <unordered_map>
#include <vector>

/**
 * @brief MPI 版 ConcatPoissonSolver2D
 *
 * 关键约定：
 * - 所有右端项与场的全量数据仅存于传入通信器的 rank 0（全局 root）。
 * - 每个 Poisson（包含 GMRES 内部调用）均通过 MPIPoissonSolver2D 在其专属子通信器中执行，
 *   子通信器的 rank 0 拥有该域的右端项，并在求解结束后将解返回给全局 root。
 * - 并行结构分层分配：root 层（tree_root）使用全部进程；其余层按该层域数量对总进程数进行近似均分，
 *   每个域获得一段连续的进程区间，构建自己的子通信器并并行求解。
 */
class MPIConcatPoissonSolver2D
{
public:
    MPIConcatPoissonSolver2D(Variable* in_variable, MPI_Comm in_comm, EnvironmentConfig* in_env_config = nullptr);
    ~MPIConcatPoissonSolver2D();

    void solve();

private:
    // 输入
    Variable*          variable = nullptr;
    EnvironmentConfig* env_config = nullptr;

    // 拓扑结构
    Geometry2D*                                                geometry = nullptr;
    Domain2DUniform*                                           tree_root = nullptr;
    std::unordered_map<Domain2DUniform*, std::unordered_map<LocationType, Domain2DUniform*>> tree_map;
    std::unordered_map<Domain2DUniform*, std::pair<LocationType, Domain2DUniform*>>           parent_map;
    std::unordered_map<Domain2DUniform*, field2*>                                      field_map;

    // 层级
    std::vector<std::vector<Domain2DUniform*>> levels;

    // MPI
    MPI_Comm global_comm = MPI_COMM_NULL;
    int      world_rank = -1, world_size = 0;

    // 每域的子通信器与求解器
    std::unordered_map<Domain2DUniform*, MPI_Comm>               domain_comm_map;
    std::unordered_map<Domain2DUniform*, DomainSolver2D*>        solver_map;
    std::unordered_map<Domain2DUniform*, field2*>                temp_fields;

private:
    // 层级构建与通信器划分
    void build_levels();
    void build_domain_communicators();
    static void even_split_ranges(int nProc, int nParts, std::vector<int>& counts, std::vector<int>& displs);

    // 构造 MPI 版 solver_map（叶子用 MPIPoisson；非叶子用 MPIGMRES）
    void construct_solver_map();

    // 根-支求解流程（仿照串行 ConcatPoissonSolver2D::solve）
    void solve_right_hand_construction();
    void solve_root_equation();
    void solve_branch_equations();

    // 数据在全局 root 与子通信器 root 之间移动
    int  comm_rank_of_local_root(MPI_Comm subcomm) const;
    void push_field_to_comm_root(MPI_Comm subcomm, field2& f_global) const;
    void pull_field_from_comm_root(MPI_Comm subcomm, field2& f_global) const;
};



#include "concat_solver2d_slab_x.h"
#include "base/domain/domain2d_mpi.h"
#include "base/parallel/mpi/distribute_slab.h"
#include "gmres_solver2d_slab_x.h"
#include "instrumentor/timer.h"
#include "io/csv_writer_2d.h"
#include "pe/poisson/poisson_solver2d_slab_x.h"

#include <string>

ConcatPoissonSolver2DSlabX::ConcatPoissonSolver2DSlabX(Variable2DSlabX* in_variable)
{
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);

    variable = in_variable;

    init_before_constructing_solver(variable);
    construct_solver_map();
}

ConcatPoissonSolver2DSlabX::~ConcatPoissonSolver2DSlabX()
{
    for (auto kv : buffer_map)
        delete[] kv.second;
}

std::unordered_map<Domain2DUniform*, DomainSolver2D*>
ConcatPoissonSolver2DSlabX::construct_local_solver_map(Domain2DUniform* domain, bool is_local, MPI_Comm local_comm)
{
    // assume domain A in layer n is adjacented with domain B, C, D in layer n+1
    // and domain A occupies rank [m, m+a], domain B occupies rank [k, k+b]
    // [m, m+a] is inrelevant with [k, k+b], they may be overlapping or not.
    // so when creating solver in A, assuming solver in B, C, D is created. we still need to create
    // solver using B, C, D domain size and A rank for A
    std::unordered_map<Domain2DUniform*, DomainSolver2D*> local_solver_map;

    // fill local_solver_map
    std::unordered_map<LocationType, Domain2DUniform*>& local_adjacency = tree_map[domain];
    for (auto kv : local_adjacency)
    {
        Domain2DUniformMPI* adjacented_domain = static_cast<Domain2DUniformMPI*>(kv.second);
        int                 adjacented_uuid   = adjacented_domain->get_uuid();

        // Construct gmres solvers in adjacented_domain but
        // communicator who belongs to current iterated domain
        if (tree_map[adjacented_domain].size() > 0)
        {
            GMRESSolver2DSlabX* adjacented_gmres_local = nullptr;
            GMRESSolver2DSlabX* adjacented_gmres       = nullptr;

            // filter process who belongs to current iterated domain
            if (is_local)
            {
                adjacented_gmres_local = new GMRESSolver2DSlabX(adjacented_domain, m, tol, maxIter, local_comm);
                adjacented_gmres_local->set_solver(new PoissonSolver2DSlabX(adjacented_domain, variable, local_comm));
            }

            // filter process who belongs to adjacented_domain
            if (variable->slab_parent_to_level.find(adjacented_domain->get_uuid()) !=
                variable->slab_parent_to_level.end())
            {
                adjacented_gmres = static_cast<GMRESSolver2DSlabX*>(solver_map[adjacented_domain]);
            }

            migrate_from(adjacented_gmres, adjacented_gmres_local);

            // filter process who belongs to current iterated domain
            if (is_local)
                local_solver_map[adjacented_domain] = adjacented_gmres_local;
        }
        // Construct pe solvers in adjacented_domain but
        // communicator who belongs to current iterated domain
        else
        {
            // filter process who belongs to current iterated domain
            if (is_local)
                local_solver_map[adjacented_domain] = new PoissonSolver2DSlabX(adjacented_domain, variable, local_comm);
        }
    }

    return local_solver_map;
}

void ConcatPoissonSolver2DSlabX::construct_solver_map_at_domain(Domain2DUniformMPI* domain)
{
    bool     is_local = false;
    MPI_Comm local_comm;
    // filter process who belongs to current iterated domain
    if (variable->slab_parent_to_level.find(domain->get_uuid()) != variable->slab_parent_to_level.end())
    {
        is_local = true;

        int local_level = variable->slab_parent_to_level[domain->get_uuid()];
        local_comm      = variable->hierarchical_slab_comms[local_level];
    }

    if (tree_map[domain].size() > 0)
    {
        GMRESSolver2DSlabX* gmres = nullptr;

        // filter process who belongs to current iterated domain
        if (is_local)
        {
            solver_map[domain] = new GMRESSolver2DSlabX(domain, m, tol, maxIter, local_comm);
            gmres              = static_cast<GMRESSolver2DSlabX*>(solver_map[domain]);
            gmres->set_solver(new PoissonSolver2DSlabX(domain, variable, local_comm));
        }

        std::unordered_map<Domain2DUniform*, DomainSolver2D*> local_solver_map =
            construct_local_solver_map(domain, is_local, local_comm);

        // filter process who belongs to current iterated domain
        if (is_local)
            gmres->schur_mat_construct(tree_map[domain], local_solver_map);

        for (auto kv : local_solver_map)
            delete kv.second;
    }
    else
    {
        // filter process who belongs to current iterated domain
        if (is_local)
            solver_map[domain] = new PoissonSolver2DSlabX(domain, variable, local_comm);
    }
}

void ConcatPoissonSolver2DSlabX::construct_solver_map()
{
    SCOPE_TIMER("[Concat] Schur complement build total = ",
                TimeRecordType::None,
                EnvironmentConfig::Get().track_pe_construct_time);

    for (auto it = hierarchical_solve_levels.rbegin(); it != hierarchical_solve_levels.rend(); ++it)
        for (Domain2DUniform* domain : *it)
            construct_solver_map_at_domain(static_cast<Domain2DUniformMPI*>(domain));

    construct_solver_map_at_domain(static_cast<Domain2DUniformMPI*>(tree_root));
}

void ConcatPoissonSolver2DSlabX::solve()
{
    EnvironmentConfig& env_cfg = EnvironmentConfig::Get();

    // Boundary
    boundary_assembly();

    // *hx*hx for each field
    for (auto kv : field_map)
    {
        auto  domain = kv.first;
        auto& f      = *kv.second;
        f *= domain->hx * domain->hx;
    }

    // Righthand construction
    int local_level = variable->hierarchical_slab_parents.size() - 1;
    for (auto it = hierarchical_solve_levels.rbegin(); it != hierarchical_solve_levels.rend(); ++it)
    {
        auto level = *it;
        for (Domain2DUniform* _domain : level)
        {
            Domain2DUniformMPI* domain = static_cast<Domain2DUniformMPI*>(_domain);

            bool is_local =
                variable->slab_parent_to_level.find(domain->get_uuid()) != variable->slab_parent_to_level.end();
            if (is_local)
                (*temp_fields[domain]) = (*field_map[domain]);

            for (auto& [location, _child_domain] : tree_map[domain])
            {
                Domain2DUniformMPI* child_domain = static_cast<Domain2DUniformMPI*>(_child_domain);

                bool is_child = variable->slab_parent_to_level.find(child_domain->get_uuid()) !=
                                variable->slab_parent_to_level.end();

                field2* f_temp       = is_local ? temp_fields[domain] : nullptr;
                field2* f            = is_local ? field_map[domain] : nullptr;
                field2* f_temp_child = is_child ? temp_fields[child_domain] : nullptr;

                bond_add_slab(child_domain, domain, location, -1.0, f_temp_child, {f_temp, f});
            }
        }

        if (local_level < 0)
        {
            std::cerr << "local_level = " << local_level << std::endl;
            continue;
        }

        Domain2DUniform* domain = variable->hierarchical_slab_parents[local_level];
        local_level--;
        {
            SCOPE_TIMER(env_cfg.pe_solve_total_name, TimeRecordType::Accumulate, false);
            solver_map[domain]->solve(*temp_fields[domain]);
        }
    }

    // Root equation
    for (auto& [location, _child_domain] : tree_map[tree_root])
    {
        Domain2DUniformMPI* child_domain = static_cast<Domain2DUniformMPI*>(_child_domain);

        bool is_child =
            variable->slab_parent_to_level.find(child_domain->get_uuid()) != variable->slab_parent_to_level.end();

        field2* f            = field_map[tree_root];
        field2* f_temp_child = is_child ? temp_fields[child_domain] : nullptr;

        bond_add_slab(child_domain, static_cast<Domain2DUniformMPI*>(tree_root), location, -1.0, f_temp_child, {f});
    }
    {
        SCOPE_TIMER(env_cfg.pe_solve_total_name, TimeRecordType::Accumulate, false);
        solver_map[tree_root]->solve(*field_map[tree_root]);
    }

    // Branch equations
    local_level = 1;
    for (auto level : hierarchical_solve_levels)
    {
        for (Domain2DUniform* _domain : level)
        {
            Domain2DUniformMPI* domain   = static_cast<Domain2DUniformMPI*>(_domain);
            LocationType        location = parent_map[domain].first;
            Domain2DUniformMPI* parent   = static_cast<Domain2DUniformMPI*>(parent_map[domain].second);

            bool is_parent =
                variable->slab_parent_to_level.find(parent->get_uuid()) != variable->slab_parent_to_level.end();

            field2* f        = field_map[domain];
            field2* f_parent = is_parent ? field_map[parent] : nullptr;

            bond_add_slab(parent, domain, location, -1.0, f_parent, {f});
        }

        Domain2DUniform* domain = variable->hierarchical_slab_parents[local_level];
        local_level++;

        {

            SCOPE_TIMER(env_cfg.pe_solve_total_name, TimeRecordType::Accumulate, false);
            solver_map[domain]->solve(*field_map[domain]);
        }
    }
}

void ConcatPoissonSolver2DSlabX::boundary_assembly()
{
    // Apply boundary conditions to all domains in the geometry
    for (auto kv : field_map)
    {
        Domain2DUniform* domain = kv.first;
        field2&          f      = *kv.second;

        if (variable->has_boundary_value_map.find(domain) == variable->has_boundary_value_map.end())
            continue;

        auto& var_has_map   = variable->has_boundary_value_map[domain];
        auto& var_value_map = variable->boundary_value_map[domain];
        auto& var_type_map  = variable->boundary_type_map[domain];

        int    nx = domain->get_nx();
        int    ny = domain->get_ny();
        double hx = domain->get_hx();
        double hy = domain->get_hy();

        PDEBoundaryType boundary_type_left  = var_type_map[LocationType::Left];
        PDEBoundaryType boundary_type_right = var_type_map[LocationType::Right];
        PDEBoundaryType boundary_type_down  = var_type_map[LocationType::Down];
        PDEBoundaryType boundary_type_up    = var_type_map[LocationType::Up];

        if (boundary_type_left == PDEBoundaryType::Dirichlet && var_has_map[LocationType::Left])
        {
            double* boundary_value = var_value_map[LocationType::Left];
            f.left_bond_add(-1.0 / hx / hx, boundary_value);
        }
        if (boundary_type_right == PDEBoundaryType::Dirichlet && var_has_map[LocationType::Right])
        {
            double* boundary_value = var_value_map[LocationType::Right];
            f.right_bond_add(-1.0 / hx / hx, boundary_value);
        }

        if (boundary_type_down == PDEBoundaryType::Dirichlet && var_has_map[LocationType::Down])
        {
            double* boundary_value = var_value_map[LocationType::Down];
            f.down_bond_add(-1.0 / hy / hy, boundary_value);
        }
        if (boundary_type_up == PDEBoundaryType::Dirichlet && var_has_map[LocationType::Up])
        {
            double* boundary_value = var_value_map[LocationType::Up];
            f.up_bond_add(-1.0 / hy / hy, boundary_value);
        }

        if (boundary_type_left == PDEBoundaryType::Neumann && var_has_map[LocationType::Left])
        {
            double* boundary_value = var_value_map[LocationType::Left];
            f.left_bond_add(1.0 / hx, boundary_value);
        }
        if (boundary_type_right == PDEBoundaryType::Neumann && var_has_map[LocationType::Right])
        {
            double* boundary_value = var_value_map[LocationType::Right];
            f.right_bond_add(-1.0 / hx, boundary_value);
        }
        if (boundary_type_down == PDEBoundaryType::Neumann && var_has_map[LocationType::Down])
        {
            double* boundary_value = var_value_map[LocationType::Down];
            f.down_bond_add(1.0 / hy, boundary_value);
        }
        if (boundary_type_up == PDEBoundaryType::Neumann && var_has_map[LocationType::Up])
        {
            double* boundary_value = var_value_map[LocationType::Up];
            f.up_bond_add(-1.0 / hy, boundary_value);
        }
    }
}

void ConcatPoissonSolver2DSlabX::bond_add_slab(Domain2DUniformMPI*         domain_src,
                                               Domain2DUniformMPI*         domain_dest,
                                               LocationType                location,
                                               double                      coeff,
                                               field2*                     f_src,
                                               const std::vector<field2*>& f_dest)
{
    bool is_dest = variable->slab_parent_to_level.find(domain_dest->get_uuid()) != variable->slab_parent_to_level.end();

    MPI_Comm comm_dest     = MPI_COMM_NULL;
    int      mpi_rank_dest = -1;
    int      mpi_size_dest = -1;
    if (is_dest)
    {
        int local_level = variable->slab_parent_to_level[domain_dest->get_uuid()];
        comm_dest       = variable->hierarchical_slab_comms[local_level];

        MPI_Comm_rank(comm_dest, &mpi_rank_dest);
        MPI_Comm_size(comm_dest, &mpi_size_dest);
    }

    bool is_src = variable->slab_parent_to_level.find(domain_src->get_uuid()) != variable->slab_parent_to_level.end();

    MPI_Comm comm_src     = MPI_COMM_NULL;
    int      mpi_rank_src = -1;
    int      mpi_size_src = -1;
    if (is_src)
    {
        int local_level_src = variable->slab_parent_to_level[domain_src->get_uuid()];
        comm_src            = variable->hierarchical_slab_comms[local_level_src];

        MPI_Comm_rank(comm_src, &mpi_rank_src);
        MPI_Comm_size(comm_src, &mpi_size_src);
    }

    if (location == LocationType::Left || location == LocationType::Right)
    {
        bool is_desired_slab_dest = false, is_desired_slab_src = false;
        if (location == LocationType::Left)
        {
            is_desired_slab_dest = mpi_rank_dest == 0;
            is_desired_slab_src  = mpi_rank_src == mpi_size_src - 1;
        }
        else
        {
            is_desired_slab_dest = mpi_rank_dest == mpi_size_dest - 1;
            is_desired_slab_src  = mpi_rank_src == 0;
        }

        int mpi_rank_src = -1, mpi_rank_dest = -1;
        // filter process who owns first slab of current iterated domain
        if (is_desired_slab_dest)
            MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank_dest);
        // filter process who owns last slab of current iterated domain's child
        if (is_desired_slab_src)
            MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank_src);

        // Synchronize metadata across MPI_COMM_WORLD
        // It should be called by all process,
        // so even though mpi_rank_src, mpi_rank_dest are only used by only one if case
        // it can not be put in if block
        MPI_Allreduce(MPI_IN_PLACE, &mpi_rank_src, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
        MPI_Allreduce(MPI_IN_PLACE, &mpi_rank_dest, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);

        // first slab of domain and last slab of child domain are located at same process
        if (is_desired_slab_dest && is_desired_slab_src)
        {
            for (auto* f : f_dest)
                f->bond_add(location, coeff, *f_src);
        }
        else
        {
            double* buffer = nullptr;
            if (is_desired_slab_dest)
                buffer = get_buffer(f_dest[0]->get_ny());

            if (is_desired_slab_src)
            {
                int i = location == LocationType::Left ? f_src->get_nx() - 1 : 0;
                MPI_Send(f_src->get_ptr(i, 0), f_src->get_ny(), MPI_DOUBLE, mpi_rank_dest, 0, MPI_COMM_WORLD);
            }
            if (is_desired_slab_dest)
                MPI_Recv(buffer, f_dest[0]->get_ny(), MPI_DOUBLE, mpi_rank_src, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

            // filter process who owns first slab of current iterated domain
            if (is_desired_slab_dest)
            {
                for (auto* f : f_dest)
                    f->bond_add(location, coeff, buffer);
            }
        }
    }
    else if (location == LocationType::Down || location == LocationType::Up)
    {
        int result = MPI_UNEQUAL;
        if (is_src && is_dest)
        {
            if (!(comm_src == MPI_COMM_NULL && comm_dest == MPI_COMM_NULL))
                MPI_Comm_compare(comm_dest, comm_src, &result);
        }

        if (result == MPI_IDENT || result == MPI_CONGRUENT)
        {
            for (auto* f : f_dest)
                f->bond_add(location, coeff, *f_src);
        }
        else
        {
            int nx_slab_dest = 0, nx_slab_src = 0;
            if (is_dest)
                nx_slab_dest = f_dest[0]->get_nx();
            if (is_src)
                nx_slab_src = f_src->get_nx();

            double *buffer_dest = nullptr, *buffer_src = nullptr;
            if (is_dest)
            {
                buffer_dest = get_buffer(nx_slab_dest);
            }
            if (is_src)
            {
                buffer_src = get_buffer(nx_slab_src);

                int j = location == LocationType::Down ? f_src->get_ny() - 1 : 0;
                for (int i = 0; i < nx_slab_src; i++)
                    buffer_src[i] = (*f_src)(i, j);
            }

            MPIUtils::redistribute_slab(buffer_src, buffer_dest, nx_slab_src, nx_slab_dest, comm_src, comm_dest);

            if (is_dest)
            {
                for (auto* f : f_dest)
                    f->bond_add(location, coeff, buffer_dest);
            }
        }
    }
}

double* ConcatPoissonSolver2DSlabX::get_buffer(int size)
{
    if (size <= 0)
        return nullptr;

    if (buffer_map.find(size) != buffer_map.end())
        return buffer_map[size];

    buffer_map[size] = new double[size];
    return buffer_map[size];
}
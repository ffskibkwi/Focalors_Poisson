#include "concat_solver2d_slab_x.h"
#include "base/domain/domain2d_mpi.h"
#include "base/parallel/mpi/distribute_slab.h"
#include "gmres_solver2d_slab_x.h"
#include "io/csv_writer_2d.h"
#include "pe/poisson/poisson_solver2d_slab_x.h"

#include <string>

ConcatPoissonSolver2DSlabX::ConcatPoissonSolver2DSlabX(Variable2DSlabX* in_variable, EnvironmentConfig* in_env_config)
{
    variable   = in_variable;
    env_config = in_env_config;

    init_before_constructing_solver(variable);
    construct_solver_map();
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
                adjacented_gmres_local =
                    new GMRESSolver2DSlabX(adjacented_domain, m, tol, maxIter, env_config, local_comm);
                adjacented_gmres_local->set_solver(
                    new PoissonSolver2DSlabX(adjacented_domain, variable, env_config, local_comm));
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
                local_solver_map[adjacented_domain] =
                    new PoissonSolver2DSlabX(adjacented_domain, variable, env_config, local_comm);
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
            solver_map[domain] = new GMRESSolver2DSlabX(domain, m, tol, maxIter, env_config, local_comm);
            gmres              = static_cast<GMRESSolver2DSlabX*>(solver_map[domain]);
            gmres->set_solver(new PoissonSolver2DSlabX(domain, variable, env_config, local_comm));
        }

        std::unordered_map<Domain2DUniform*, DomainSolver2D*> local_solver_map =
            construct_local_solver_map(domain, is_local, local_comm);

        // filter process who belongs to current iterated domain
        if (is_local)
        {
            std::chrono::steady_clock::time_point t0, t1;
            if (track_time)
                t0 = std::chrono::steady_clock::now();

            gmres->schur_mat_construct(tree_map[domain], local_solver_map);

            if (track_time)
            {
                t1 = std::chrono::steady_clock::now();
                schur_total += t1 - t0;
            }
        }
    }
    else
    {
        // filter process who belongs to current iterated domain
        if (is_local)
            solver_map[domain] = new PoissonSolver2DSlabX(domain, variable, env_config, local_comm);
    }
}

void ConcatPoissonSolver2DSlabX::construct_solver_map()
{
    if (track_time)
        schur_total = std::chrono::duration<double, std::milli>(0.0);

    for (auto it = hierarchical_solve_levels.rbegin(); it != hierarchical_solve_levels.rend(); ++it)
        for (Domain2DUniform* domain : *it)
            construct_solver_map_at_domain(static_cast<Domain2DUniformMPI*>(domain));

    construct_solver_map_at_domain(static_cast<Domain2DUniformMPI*>(tree_root));

    if (track_time)
        std::cout << "[Concat] Schur complement build total=" << schur_total.count() << " ms" << std::endl;
}

void ConcatPoissonSolver2DSlabX::solve()
{
    const bool track_detail_time = env_config && env_config->track_pe_solve_detail_time;
    const bool track_total_time  = env_config && env_config->track_pe_solve_total_time;

    // Boundary
    boundary_assembly();

    // *hx*hx for each field
    for (auto kv : field_map)
    {
        auto  domain = kv.first;
        auto& f      = *kv.second;
        f            = f * (domain->hx * domain->hx);
    }

    auto solve_start = std::chrono::steady_clock::time_point();
    if (track_total_time)
        solve_start = std::chrono::steady_clock::now();

    // Righthand construction
    for (auto it = hierarchical_solve_levels.rbegin(); it != hierarchical_solve_levels.rend(); ++it)
    {
        auto level = *it;
        for (Domain2DUniform* _domain : level)
        {
            Domain2DUniformMPI* domain = static_cast<Domain2DUniformMPI*>(_domain);

            bool is_local =
                variable->slab_parent_to_level.find(domain->get_uuid()) != variable->slab_parent_to_level.end();

            MPI_Comm comm_local;
            int      mpi_rank_local = -1;
            int      mpi_size_local = -1;
            if (is_local)
            {
                int local_level = variable->slab_parent_to_level[domain->get_uuid()];
                comm_local      = variable->hierarchical_slab_comms[local_level];

                MPI_Comm_rank(comm_local, &mpi_rank_local);
                MPI_Comm_size(comm_local, &mpi_size_local);
            }

            // filter process who belongs to current iterated domain
            if (is_local)
                (*temp_fields[domain]) = (*field_map[domain]);

            for (auto& [location, _child_domain] : tree_map[domain])
            {
                Domain2DUniformMPI* child_domain = static_cast<Domain2DUniformMPI*>(_child_domain);

                bool is_child = variable->slab_parent_to_level.find(child_domain->get_uuid()) !=
                                variable->slab_parent_to_level.end();

                MPI_Comm comm_child;
                int      mpi_rank_child = -1;
                int      mpi_size_child = -1;
                // filter process who belongs to current iterated domain's child
                if (is_child)
                {
                    int local_level_child = variable->slab_parent_to_level[child_domain->get_uuid()];
                    comm_child            = variable->hierarchical_slab_comms[local_level_child];

                    MPI_Comm_rank(comm_child, &mpi_rank_child);
                    MPI_Comm_size(comm_child, &mpi_size_child);
                }

                if (location == LocationType::Left || location == LocationType::Right)
                {
                    bool is_desired_slab_local = false, is_desired_slab_child = false;
                    if (location == LocationType::Left)
                    {
                        is_desired_slab_local = mpi_rank_local == 0;
                        is_desired_slab_child = mpi_rank_child == mpi_size_child - 1;
                    }
                    else
                    {
                        is_desired_slab_local = mpi_rank_local == mpi_size_local - 1;
                        is_desired_slab_child = mpi_rank_child == 0;
                    }

                    // first slab of domain and last slab of child domain are located at same process
                    if (is_desired_slab_local && is_desired_slab_child)
                    {
                        temp_fields[domain]->bond_add(location, -1., *temp_fields[child_domain]);
                        field_map[domain]->bond_add(location, -1., *temp_fields[child_domain]);
                    }
                    else
                    {
                        double* buffer = nullptr;
                        if (is_desired_slab_local)
                            buffer = get_buffer(temp_fields[domain]->get_ny());

                        int mpi_rank_src, mpi_rank_dest;
                        // filter process who owns first slab of current iterated domain
                        if (is_desired_slab_local)
                            MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank_dest);
                        // filter process who owns last slab of current iterated domain's child
                        if (is_desired_slab_child)
                            MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank_src);

                        // Synchronize metadata across MPI_COMM_WORLD
                        MPI_Allreduce(MPI_IN_PLACE, &mpi_rank_src, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
                        MPI_Allreduce(MPI_IN_PLACE, &mpi_rank_dest, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);

                        if (is_desired_slab_child)
                        {
                            int i = location == LocationType::Left ? temp_fields[child_domain]->get_nx() - 1 : 0;
                            MPI_Send(temp_fields[child_domain]->get_ptr(i, 0),
                                     temp_fields[child_domain]->get_ny(),
                                     MPI_DOUBLE,
                                     mpi_rank_dest,
                                     0,
                                     MPI_COMM_WORLD);
                        }
                        if (is_desired_slab_local)
                            MPI_Recv(buffer,
                                     temp_fields[domain]->get_ny(),
                                     MPI_DOUBLE,
                                     mpi_rank_src,
                                     0,
                                     MPI_COMM_WORLD,
                                     MPI_STATUS_IGNORE);

                        // filter process who owns first slab of current iterated domain
                        if (is_desired_slab_local)
                        {
                            temp_fields[domain]->bond_add(location, -1., buffer);
                            field_map[domain]->bond_add(location, -1., buffer);
                        }
                    }
                }
                else if (location == LocationType::Down || location == LocationType::Up)
                {
                    int result;
                    MPI_Comm_compare(comm_local, comm_child, &result);
                    if (result == MPI_IDENT || result == MPI_CONGRUENT)
                    {
                        temp_fields[domain]->bond_add(location, -1., *temp_fields[child_domain]);
                        field_map[domain]->bond_add(location, -1., *temp_fields[child_domain]);
                    }
                    else
                    {
                        int nx_slab_local = 0, nx_slab_child = 0;
                        if (is_local)
                            nx_slab_local = temp_fields[domain]->get_nx();
                        if (is_child)
                            nx_slab_child = temp_fields[child_domain]->get_nx();

                        double *buffer_local = nullptr, *buffer_child = nullptr;
                        if (is_local)
                        {
                            buffer_local = get_buffer(nx_slab_local);
                        }
                        if (is_child)
                        {
                            buffer_child = get_buffer(nx_slab_child);

                            int j = location == LocationType::Down ? 0 : temp_fields[child_domain]->get_ny() - 1;
                            for (int i = 0; i < nx_slab_child; i++)
                                buffer_child[i] = (*temp_fields[child_domain])(i, j);
                        }

                        MPIUtils::redistribute_slab(
                            buffer_child, buffer_local, nx_slab_child, nx_slab_local, comm_child, comm_local);

                        if (is_local)
                        {
                            temp_fields[domain]->bond_add(location, -1., buffer_local);
                            field_map[domain]->bond_add(location, -1., buffer_local);
                        }
                    }
                }
            }

            if (env_config && env_config->showCurrentStep)
            {
                double s_pre = temp_fields[domain]->sum();
                std::cout << "[Concat] Domain " << domain->name << " temp sum before solve=" << s_pre << std::endl;
            }

            if (track_detail_time)
            {
                auto t0 = std::chrono::steady_clock::now();
                solver_map[domain]->solve(*temp_fields[domain]);
                auto t1 = std::chrono::steady_clock::now();
                std::cout << "[Concat] Domain " << domain->name
                          << " solve time=" << std::chrono::duration<double, std::milli>(t1 - t0).count() << " ms"
                          << std::endl;
            }
            else
            {
                solver_map[domain]->solve(*temp_fields[domain]);
            }

            if (env_config && env_config->debug_concat)
            {
                std::string fname_Ainv = env_config->debugOutputDir + "/Ainv_f_" + domain->name;
                IO::write_csv(*temp_fields[domain], fname_Ainv);
            }
            if (env_config && env_config->showCurrentStep)
            {
                double s_post = temp_fields[domain]->sum();
                std::cout << "[Concat] Domain " << domain->name << " temp sum after solve=" << s_post << std::endl;
            }
        }
    }

    // Root equation
    for (auto& [location, child_domain] : tree_map[tree_root])
        field_map[tree_root]->bond_add(location, -1., *temp_fields[child_domain]);

    if (env_config && env_config->showCurrentStep)
    {
        double s_root_pre = field_map[tree_root]->sum();
        std::cout << "[Concat] Root domain " << tree_root->name << " field sum before solve=" << s_root_pre
                  << std::endl;
    }

    if (track_detail_time)
    {
        auto t0 = std::chrono::steady_clock::now();
        solver_map[tree_root]->solve(*field_map[tree_root]);
        auto t1 = std::chrono::steady_clock::now();
        std::cout << "[Concat] Root domain " << tree_root->name
                  << " solve time=" << std::chrono::duration<double, std::milli>(t1 - t0).count() << " ms" << std::endl;
    }
    else
    {
        solver_map[tree_root]->solve(*field_map[tree_root]);
    }

    if (env_config && env_config->showCurrentStep)
    {
        double s_root_post = field_map[tree_root]->sum();
        std::cout << "[Concat] Root domain " << tree_root->name << " field sum after solve=" << s_root_post
                  << std::endl;
    }

    // Branch equations
    for (auto level : hierarchical_solve_levels)
    {
        for (Domain2DUniform* domain : level)
        {
            field_map[domain]->bond_add(parent_map[domain].first, -1., *field_map[parent_map[domain].second]);

            if (env_config && env_config->showCurrentStep)
            {
                double s_pre = field_map[domain]->sum();
                std::cout << "[Concat] Branch domain " << domain->name << " field sum before solve=" << s_pre
                          << std::endl;
            }

            if (track_detail_time)
            {
                auto t0 = std::chrono::steady_clock::now();
                solver_map[domain]->solve(*field_map[domain]);
                auto t1 = std::chrono::steady_clock::now();
                std::cout << "[Concat] Branch domain " << domain->name
                          << " solve time=" << std::chrono::duration<double, std::milli>(t1 - t0).count() << " ms"
                          << std::endl;
            }
            else
            {
                solver_map[domain]->solve(*field_map[domain]);
            }

            if (env_config && env_config->showCurrentStep)
            {
                double s_post = field_map[domain]->sum();
                std::cout << "[Concat] Branch domain " << domain->name << " field sum after solve=" << s_post
                          << std::endl;
            }
        }
    }

    if (track_total_time)
    {
        auto solve_end = std::chrono::steady_clock::now();
        std::cout << "[Concat] Solve total (exclude boundary/scale)="
                  << std::chrono::duration<double, std::milli>(solve_end - solve_start).count() << " ms" << std::endl;
    }
}

void ConcatPoissonSolver2DSlabX::boundary_assembly()
{
    // Apply boundary conditions to all domains in the geometry
    for (auto& domain : variable->geometry->domains)
    {
        if (variable->has_boundary_value_map.find(domain) == variable->has_boundary_value_map.end())
            continue;

        auto& var_has_map   = variable->has_boundary_value_map[domain];
        auto& var_value_map = variable->boundary_value_map[domain];
        auto& var_type_map  = variable->boundary_type_map[domain];

        field2& f = *field_map[domain];

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

double* ConcatPoissonSolver2DSlabX::get_buffer(int size)
{
    if (size <= 0)
        return nullptr;

    if (buffer_map.find(size) != buffer_map.end())
        return buffer_map[size];

    buffer_map[size] = new double[size];
    return buffer_map[size];
}
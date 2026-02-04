#include "concat_solver2d.h"
#include "gmres_solver2d.h"
#include "io/csv_writer_2d.h"
#include "pe/poisson/poisson_solver2d.h"

#include <string>

void ConcatPoissonSolver2D::set_parameter(int in_m, double in_tol, int in_maxIter)
{
    m       = in_m;
    tol     = in_tol;
    maxIter = in_maxIter;
}

ConcatPoissonSolver2D::ConcatPoissonSolver2D(Variable2D* in_variable, EnvironmentConfig* in_env_config)
    : variable(in_variable)
    , env_config(in_env_config)
{
    init_before_constructing_solver(variable);
    construct_solver_map();
}

ConcatPoissonSolver2D::~ConcatPoissonSolver2D()
{
    for (auto& [domain, temp_field] : temp_fields)
        delete temp_field;
    for (auto kv : solver_map)
        delete kv.second;
}

void ConcatPoissonSolver2D::init_before_constructing_solver(Variable2D* _variable)
{
    // geometry double check
    if (_variable->geometry == nullptr)
        throw std::runtime_error("ConcatPoissonSolver2D: _variable->geometry is null");
    if (!_variable->geometry->is_checked)
        _variable->geometry->check();
    if (_variable->geometry->tree_root == nullptr || _variable->geometry->tree_map.empty())
        _variable->geometry->solve_prepare();

    tree_root                 = _variable->geometry->tree_root;
    tree_map                  = _variable->geometry->tree_map;
    parent_map                = _variable->geometry->parent_map;
    field_map                 = _variable->field_map;
    hierarchical_solve_levels = _variable->geometry->hierarchical_solve_levels;
    hierarchical_solve_levels.erase(hierarchical_solve_levels.begin()); // erase tree_root

    // Construct the temp field for each domain
    for (auto& [domain, field] : field_map)
    {
        if (domain != tree_root)
            temp_fields[domain] = new field2(field->get_nx(), field->get_ny(), field->get_name() + "_temp");
    }

    showGmresRes = env_config && env_config->showGmresRes;
    track_time   = env_config && env_config->track_pe_construct_time;
}

void ConcatPoissonSolver2D::construct_solver_map_at_domain(Domain2DUniform* domain)
{
    if (tree_map[domain].size() > 0)
    {
        solver_map[domain] = new GMRESSolver2D(domain, m, tol, maxIter, env_config);
        auto* gmres        = static_cast<GMRESSolver2D*>(solver_map[domain]);
        gmres->set_solver(new PoissonSolver2D(domain, variable, env_config));

        std::chrono::steady_clock::time_point t0, t1;
        if (track_time)
            t0 = std::chrono::steady_clock::now();

        gmres->schur_mat_construct(tree_map[domain], solver_map);

        if (track_time)
        {
            t1 = std::chrono::steady_clock::now();
            schur_total += t1 - t0;
        }
    }
    else
    {
        solver_map[domain] = new PoissonSolver2D(domain, variable, env_config);
    }
}

void ConcatPoissonSolver2D::construct_solver_map()
{
    if (track_time)
        schur_total = std::chrono::duration<double, std::milli>(0.0);

    for (auto it = hierarchical_solve_levels.rbegin(); it != hierarchical_solve_levels.rend(); ++it)
        for (Domain2DUniform* domain : *it)
            construct_solver_map_at_domain(domain);

    construct_solver_map_at_domain(tree_root);

    if (track_time)
        std::cout << "[Concat] Schur complement build total=" << schur_total.count() << " ms" << std::endl;
}

void ConcatPoissonSolver2D::solve()
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
        for (Domain2DUniform* domain : level)
        {
            (*temp_fields[domain]) = (*field_map[domain]);
            for (auto& [location, child_domain] : tree_map[domain])
            {
                temp_fields[domain]->bond_add(location, -1., *temp_fields[child_domain]);
                field_map[domain]->bond_add(location, -1., *temp_fields[child_domain]);
            }
        }
        for (Domain2DUniform* domain : level)
        {
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
        }
        for (Domain2DUniform* domain : level)
        {
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

void ConcatPoissonSolver2D::boundary_assembly()
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

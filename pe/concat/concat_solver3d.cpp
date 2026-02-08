#include "concat_solver3d.h"
#include "gmres_solver3d.h"
#include "pe/poisson/poisson_solver3d.h"

void ConcatPoissonSolver3D::set_parameter(int in_m, double in_tol, int in_maxIter)
{
    m       = in_m;
    tol     = in_tol;
    maxIter = in_maxIter;
}

ConcatPoissonSolver3D::ConcatPoissonSolver3D(Variable3D* in_variable)
    : variable(in_variable)
{
    init_before_constructing_solver(variable);
    construct_solver_map();
}

ConcatPoissonSolver3D::~ConcatPoissonSolver3D()
{
    for (auto& [domain, temp_field] : temp_fields)
        delete temp_field;
    for (auto kv : solver_map)
        delete kv.second;
}

void ConcatPoissonSolver3D::init_before_constructing_solver(Variable3D* _variable)
{
    // geometry double check
    if (_variable->geometry == nullptr)
        throw std::runtime_error("ConcatPoissonSolver3D: _variable->geometry is null");
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
            temp_fields[domain] =
                new field3(field->get_nx(), field->get_ny(), field->get_nz(), field->get_name() + "_temp");
    }

    EnvironmentConfig& env_cfg = EnvironmentConfig::Get();

    showGmresRes = env_cfg.showGmresRes;
    track_time   = env_cfg.track_pe_construct_time;
}

void ConcatPoissonSolver3D::construct_solver_map_at_domain(Domain3DUniform* domain)
{
    if (tree_map[domain].size() > 0)
    {
        solver_map[domain] = new GMRESSolver3D(domain, m, tol, maxIter);
        auto* gmres        = static_cast<GMRESSolver3D*>(solver_map[domain]);
        gmres->schur_mat_construct(tree_map[domain], solver_map);
        gmres->set_solver(new PoissonSolver3D(domain, variable));
    }
    else
    {
        solver_map[domain] = new PoissonSolver3D(domain, variable);
    }
}

void ConcatPoissonSolver3D::construct_solver_map()
{
    for (auto it = hierarchical_solve_levels.rbegin(); it != hierarchical_solve_levels.rend(); ++it)
        for (Domain3DUniform* domain : *it)
            construct_solver_map_at_domain(domain);

    construct_solver_map_at_domain(tree_root);
}

void ConcatPoissonSolver3D::solve()
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

    // Righthand construction (bottom-up pass)
    if (env_cfg.showCurrentStep)
        std::cout << "[Concat3D] Bottom-up pass: start" << std::endl;
    for (auto it = hierarchical_solve_levels.rbegin(); it != hierarchical_solve_levels.rend(); ++it)
    {
        auto level = *it;
        for (Domain3DUniform* domain : level)
        {
            if (env_cfg.showCurrentStep)
                std::cout << "[Concat3D] Bottom-up pass: solve domain " << domain->name << std::endl;
            (*temp_fields[domain]) = (*field_map[domain]);
            for (auto& [location, child_domain] : tree_map[domain])
            {
                temp_fields[domain]->bond_add(location, -1., *temp_fields[child_domain]);
                field_map[domain]->bond_add(location, -1., *temp_fields[child_domain]);
            }
        }
        for (Domain3DUniform* domain : level)
        {
            if (env_cfg.showCurrentStep)
            {
                double s_pre = temp_fields[domain]->sum();
                std::cout << "[Concat3D] Bottom-up pass: domain " << domain->name << " temp sum before solve=" << s_pre
                          << std::endl;
            }

            solver_map[domain]->solve(*temp_fields[domain]);

            if (env_cfg.showCurrentStep)
            {
                double s_post = temp_fields[domain]->sum();
                std::cout << "[Concat3D] Bottom-up pass: domain " << domain->name << " temp sum after solve=" << s_post
                          << std::endl;
            }
        }
    }
    if (env_cfg.showCurrentStep)
        std::cout << "[Concat3D] Bottom-up pass: done" << std::endl;

    // Root equation
    for (auto& [location, child_domain] : tree_map[tree_root])
        field_map[tree_root]->bond_add(location, -1., *temp_fields[child_domain]);

    if (env_cfg.showCurrentStep)
    {
        double s_root_pre = field_map[tree_root]->sum();
        std::cout << "[Concat3D] Root solve: domain " << tree_root->name << " field sum before solve=" << s_root_pre
                  << std::endl;
    }

    solver_map[tree_root]->solve(*field_map[tree_root]);

    if (env_cfg.showCurrentStep)
    {
        double s_root_post = field_map[tree_root]->sum();
        std::cout << "[Concat3D] Root solve: domain " << tree_root->name << " field sum after solve=" << s_root_post
                  << std::endl;
    }

    // Branch equations (top-down pass)
    if (env_cfg.showCurrentStep)
        std::cout << "[Concat3D] Top-down pass: start" << std::endl;
    for (auto level : hierarchical_solve_levels)
    {
        for (Domain3DUniform* domain : level)
        {
            if (env_cfg.showCurrentStep)
            {
                auto parent = parent_map[domain];
                std::cout << "[Concat3D] Top-down pass: solve domain " << domain->name << " (parent "
                          << parent.second->name << " at " << parent.first << ")" << std::endl;
            }

            field_map[domain]->bond_add(parent_map[domain].first, -1., *field_map[parent_map[domain].second]);
        }
        for (Domain3DUniform* domain : level)
        {
            if (env_cfg.showCurrentStep)
            {
                double s_pre = field_map[domain]->sum();
                std::cout << "[Concat3D] Top-down pass: domain " << domain->name << " field sum before solve=" << s_pre
                          << std::endl;
            }

            solver_map[domain]->solve(*field_map[domain]);

            if (env_cfg.showCurrentStep)
            {
                double s_post = field_map[domain]->sum();
                std::cout << "[Concat3D] Top-down pass: domain " << domain->name << " field sum after solve=" << s_post
                          << std::endl;
            }
        }
    }
    if (env_cfg.showCurrentStep)
        std::cout << "[Concat3D] Top-down pass: done" << std::endl;
}

void ConcatPoissonSolver3D::boundary_assembly()
{
    // Apply boundary conditions to all domains in the geometry
    for (auto& domain : variable->geometry->domains)
    {
        if (variable->has_boundary_value_map.find(domain) == variable->has_boundary_value_map.end())
            continue;

        auto& var_has_map   = variable->has_boundary_value_map[domain];
        auto& var_value_map = variable->boundary_value_map[domain];
        auto& var_type_map  = variable->boundary_type_map[domain];

        field3& f = *field_map[domain];

        int    nx = domain->get_nx();
        int    ny = domain->get_ny();
        int    nz = domain->get_nz();
        double hx = domain->get_hx();
        double hy = domain->get_hy();
        double hz = domain->get_hz();

        PDEBoundaryType boundary_type_left  = var_type_map[LocationType::Left];
        PDEBoundaryType boundary_type_right = var_type_map[LocationType::Right];
        PDEBoundaryType boundary_type_front = var_type_map[LocationType::Front];
        PDEBoundaryType boundary_type_back  = var_type_map[LocationType::Back];
        PDEBoundaryType boundary_type_down  = var_type_map[LocationType::Down];
        PDEBoundaryType boundary_type_up    = var_type_map[LocationType::Up];

        if (boundary_type_left == PDEBoundaryType::Dirichlet && var_has_map[LocationType::Left])
        {
            field2& boundary_value = *var_value_map[LocationType::Left];
            f.left_bond_add(-1.0 / hx / hx, boundary_value);
        }
        if (boundary_type_right == PDEBoundaryType::Dirichlet && var_has_map[LocationType::Right])
        {
            field2& boundary_value = *var_value_map[LocationType::Right];
            f.right_bond_add(-1.0 / hx / hx, boundary_value);
        }

        if (boundary_type_front == PDEBoundaryType::Dirichlet && var_has_map[LocationType::Front])
        {
            field2& boundary_value = *var_value_map[LocationType::Front];
            f.front_bond_add(-1.0 / hy / hy, boundary_value);
        }
        if (boundary_type_back == PDEBoundaryType::Dirichlet && var_has_map[LocationType::Back])
        {
            field2& boundary_value = *var_value_map[LocationType::Back];
            f.back_bond_add(-1.0 / hy / hy, boundary_value);
        }

        if (boundary_type_down == PDEBoundaryType::Dirichlet && var_has_map[LocationType::Down])
        {
            field2& boundary_value = *var_value_map[LocationType::Down];
            f.down_bond_add(-1.0 / hz / hz, boundary_value);
        }
        if (boundary_type_up == PDEBoundaryType::Dirichlet && var_has_map[LocationType::Up])
        {
            field2& boundary_value = *var_value_map[LocationType::Up];
            f.up_bond_add(-1.0 / hz / hz, boundary_value);
        }

        if (boundary_type_left == PDEBoundaryType::Neumann && var_has_map[LocationType::Left])
        {
            field2& boundary_value = *var_value_map[LocationType::Left];
            f.left_bond_add(1.0 / hx, boundary_value);
        }
        if (boundary_type_right == PDEBoundaryType::Neumann && var_has_map[LocationType::Right])
        {
            field2& boundary_value = *var_value_map[LocationType::Right];
            f.right_bond_add(-1.0 / hx, boundary_value);
        }
        if (boundary_type_front == PDEBoundaryType::Neumann && var_has_map[LocationType::Front])
        {
            field2& boundary_value = *var_value_map[LocationType::Front];
            f.front_bond_add(1.0 / hy, boundary_value);
        }
        if (boundary_type_back == PDEBoundaryType::Neumann && var_has_map[LocationType::Back])
        {
            field2& boundary_value = *var_value_map[LocationType::Back];
            f.back_bond_add(-1.0 / hy, boundary_value);
        }
        if (boundary_type_down == PDEBoundaryType::Neumann && var_has_map[LocationType::Down])
        {
            field2& boundary_value = *var_value_map[LocationType::Down];
            f.down_bond_add(1.0 / hz, boundary_value);
        }
        if (boundary_type_up == PDEBoundaryType::Neumann && var_has_map[LocationType::Up])
        {
            field2& boundary_value = *var_value_map[LocationType::Up];
            f.up_bond_add(-1.0 / hz, boundary_value);
        }
    }
}

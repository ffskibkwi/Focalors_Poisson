#include "concat_solver3d.h"
#include "base/domain/geometry_tree.hpp"

namespace
{
    const char* location_to_string(LocationType location)
    {
        switch (location)
        {
            case LocationType::Left:
                return "Left";
            case LocationType::Right:
                return "Right";
            case LocationType::Front:
                return "Front";
            case LocationType::Back:
                return "Back";
            case LocationType::Down:
                return "Down";
            case LocationType::Up:
                return "Up";
            default:
                return "Unknown";
        }
    }
} // namespace

void ConcatPoissonSolver3D::set_parameter(int in_m, double in_tol, int in_maxIter)
{
    m       = in_m;
    tol     = in_tol;
    maxIter = in_maxIter;
}

ConcatPoissonSolver3D::ConcatPoissonSolver3D(Variable3D* in_variable, EnvironmentConfig* in_env_config)
    : variable(in_variable)
    , env_config(in_env_config)
{
    // config load
    if (in_env_config)
    {
        showGmresRes = in_env_config->showGmresRes;
    }

    // geometry double check
    if (variable->geometry == nullptr)
        throw std::runtime_error("ConcatPoissonSolver3D: variable->geometry is null");
    if (!variable->geometry->is_checked)
        variable->geometry->check();
    if (variable->geometry->tree_root == nullptr || variable->geometry->tree_map.empty())
        variable->geometry->solve_prepare();

    tree_root  = variable->geometry->tree_root;
    tree_map   = variable->geometry->tree_map;
    parent_map = variable->geometry->parent_map;
    field_map  = variable->field_map;

    solve_order = TreeUtils::buildLevelsFromTree(tree_root, tree_map);
    solve_order.erase(solve_order.begin()); // erase tree_root
    construct_solver_map();

    // Construct the temp field for each domain
    for (auto& [domain, field] : field_map)
    {
        if (domain != tree_root)
            temp_fields[domain] =
                new field3(field->get_nx(), field->get_ny(), field->get_nz(), field->get_name() + "_temp");
    }
}

ConcatPoissonSolver3D::~ConcatPoissonSolver3D()
{
    for (auto& [domain, temp_field] : temp_fields)
        delete temp_field;
    for (auto kv : solver_map)
        delete kv.second;
}

void ConcatPoissonSolver3D::construct_solver_map()
{
    // Construct solvers (for non-root domains)
    for (auto it = solve_order.rbegin(); it != solve_order.rend(); ++it)
    {
        auto level = *it;
        for (Domain3DUniform* domain : level)
        {
            if (tree_map[domain].size() > 0)
            {
                solver_map[domain] = new GMRESSolver3D(domain, variable, m, tol, maxIter, env_config);
                static_cast<GMRESSolver3D*>(solver_map[domain])->schur_mat_construct(tree_map[domain], solver_map);
            }
            else
            {
                solver_map[domain] = new PoissonSolver3D(domain, variable, env_config);
            }
        }
    }

    // Construct solver for root domain
    if (tree_map[tree_root].size() > 0)
    {
        solver_map[tree_root] = new GMRESSolver3D(tree_root, variable, m, tol, maxIter, env_config);
        static_cast<GMRESSolver3D*>(solver_map[tree_root])->schur_mat_construct(tree_map[tree_root], solver_map);
    }
    else
    {
        solver_map[tree_root] = new PoissonSolver3D(tree_root, variable, env_config);
    }
}

void ConcatPoissonSolver3D::solve()
{
    // Boundary
    boundary_assembly();

    // *hx*hx for each field
    for (auto& domain : variable->geometry->domains)
    {
        field3& f = *field_map[domain];
        f         = f * (domain->hx * domain->hx);
    }

    // Righthand construction (bottom-up pass)
    if (env_config && env_config->showCurrentStep)
        std::cout << "[Concat3D] Bottom-up pass: start" << std::endl;
    for (auto it = solve_order.rbegin(); it != solve_order.rend(); ++it)
    {
        auto level = *it;
        for (Domain3DUniform* domain : level)
        {
            if (env_config && env_config->showCurrentStep)
                std::cout << "[Concat3D] Bottom-up pass: solve domain " << domain->name << std::endl;
            (*temp_fields[domain]) = (*field_map[domain]);
            for (auto& [location, child_domain] : tree_map[domain])
            {
                temp_fields[domain]->bond_add(location, -1., *temp_fields[child_domain]);
                field_map[domain]->bond_add(location, -1., *temp_fields[child_domain]);
            }

            if (env_config && env_config->showCurrentStep)
            {
                double s_pre = temp_fields[domain]->sum();
                std::cout << "[Concat3D] Bottom-up pass: domain " << domain->name << " temp sum before solve=" << s_pre
                          << std::endl;
            }

            solver_map[domain]->solve(*temp_fields[domain]);

            if (env_config && env_config->showCurrentStep)
            {
                double s_post = temp_fields[domain]->sum();
                std::cout << "[Concat3D] Bottom-up pass: domain " << domain->name << " temp sum after solve=" << s_post
                          << std::endl;
            }
        }
    }
    if (env_config && env_config->showCurrentStep)
        std::cout << "[Concat3D] Bottom-up pass: done" << std::endl;

    // Root equation
    for (auto& [location, child_domain] : tree_map[tree_root])
        field_map[tree_root]->bond_add(location, -1., *temp_fields[child_domain]);

    if (env_config && env_config->showCurrentStep)
    {
        double s_root_pre = field_map[tree_root]->sum();
        std::cout << "[Concat3D] Root solve: domain " << tree_root->name << " field sum before solve=" << s_root_pre
                  << std::endl;
    }

    solver_map[tree_root]->solve(*field_map[tree_root]);

    if (env_config && env_config->showCurrentStep)
    {
        double s_root_post = field_map[tree_root]->sum();
        std::cout << "[Concat3D] Root solve: domain " << tree_root->name << " field sum after solve=" << s_root_post
                  << std::endl;
    }

    // Branch equations (top-down pass)
    if (env_config && env_config->showCurrentStep)
        std::cout << "[Concat3D] Top-down pass: start" << std::endl;
    for (auto level : solve_order)
    {
        for (Domain3DUniform* domain : level)
        {
            if (env_config && env_config->showCurrentStep)
            {
                auto parent = parent_map[domain];
                std::cout << "[Concat3D] Top-down pass: solve domain " << domain->name << " (parent "
                          << parent.second->name << " at " << location_to_string(parent.first) << ")" << std::endl;
            }

            field_map[domain]->bond_add(parent_map[domain].first, -1., *field_map[parent_map[domain].second]);

            if (env_config && env_config->showCurrentStep)
            {
                double s_pre = field_map[domain]->sum();
                std::cout << "[Concat3D] Top-down pass: domain " << domain->name << " field sum before solve=" << s_pre
                          << std::endl;
            }

            solver_map[domain]->solve(*field_map[domain]);

            if (env_config && env_config->showCurrentStep)
            {
                double s_post = field_map[domain]->sum();
                std::cout << "[Concat3D] Top-down pass: domain " << domain->name << " field sum after solve=" << s_post
                          << std::endl;
            }
        }
    }
    if (env_config && env_config->showCurrentStep)
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

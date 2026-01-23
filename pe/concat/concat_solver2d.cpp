#include "concat_solver2d.h"
#include "base/domain/geometry_tree.hpp"
#include "io/csv_writer_2d.h"

#include <string>

void ConcatPoissonSolver2D::set_parameter(int in_m, double in_tol, int in_maxIter)
{
    m       = in_m;
    tol     = in_tol;
    maxIter = in_maxIter;
}

ConcatPoissonSolver2D::ConcatPoissonSolver2D(Variable* in_variable, EnvironmentConfig* in_env_config)
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
        throw std::runtime_error("ConcatPoissonSolver2D: variable->geometry is null");
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
            temp_fields[domain] = new field2(field->get_nx(), field->get_ny(), field->get_name() + "_temp");
    }
}

ConcatPoissonSolver2D::~ConcatPoissonSolver2D()
{
    for (auto& [domain, temp_field] : temp_fields)
        delete temp_field;
    for (auto level : solve_order)
        for (Domain2DUniform* domain : level)
            delete solver_map[domain];
}

void ConcatPoissonSolver2D::construct_solver_map()
{
    // Construct solvers (for non-root domains)
    for (auto it = solve_order.rbegin(); it != solve_order.rend(); ++it)
    {
        auto level = *it;
        for (Domain2DUniform* domain : level)
        {
            if (tree_map[domain].size() > 0)
            {
                solver_map[domain] = new GMRESSolver2D(domain, variable, m, tol, maxIter, env_config);
                static_cast<GMRESSolver2D*>(solver_map[domain])
                    ->schur_mat_construct(tree_map[domain], solver_map); // Here use RTTI
            }
            else
            {
                solver_map[domain] = new PoissonSolver2D(domain, variable, env_config);
            }
        }
    }

    // Construct solver for root domain
    if (tree_map[tree_root].size() > 0)
    {
        solver_map[tree_root] = new GMRESSolver2D(tree_root, variable, m, tol, maxIter, env_config);
        static_cast<GMRESSolver2D*>(solver_map[tree_root])
            ->schur_mat_construct(tree_map[tree_root], solver_map); // Here use RTTI
    }
    else
    {
        solver_map[tree_root] = new PoissonSolver2D(tree_root, variable, env_config);
    }
}

void ConcatPoissonSolver2D::solve()
{
    // Boundary
    boundary_assembly();

    // *hx*hx for each field
    for (auto& domain : variable->geometry->domains)
    {
        field2& f = *field_map[domain];
        f         = f * (domain->hx * domain->hx);
    }

    // Righthand construction
    for (auto it = solve_order.rbegin(); it != solve_order.rend(); ++it)
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

            if (env_config && env_config->showCurrentStep)
            {
                double s_pre = temp_fields[domain]->sum();
                std::cout << "[Concat] Domain " << domain->name << " temp sum before solve=" << s_pre << std::endl;
            }

            solver_map[domain]->solve(*temp_fields[domain]);

            if (env_config && env_config->debugMode)
            {
                std::string fname_Ainv = env_config->debugOutputDir + "/Ainv_f_" + domain->name;
                IO::field_to_csv(*temp_fields[domain], fname_Ainv);
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

    solver_map[tree_root]->solve(*field_map[tree_root]);

    if (env_config && env_config->showCurrentStep)
    {
        double s_root_post = field_map[tree_root]->sum();
        std::cout << "[Concat] Root domain " << tree_root->name << " field sum after solve=" << s_root_post
                  << std::endl;
    }

    // Branch equations
    for (auto level : solve_order)
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

            solver_map[domain]->solve(*field_map[domain]);

            if (env_config && env_config->showCurrentStep)
            {
                double s_post = field_map[domain]->sum();
                std::cout << "[Concat] Branch domain " << domain->name << " field sum after solve=" << s_post
                          << std::endl;
            }
        }
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
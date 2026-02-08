#include "concat_solver2d.h"
#include "gmres_solver2d.h"
#include "instrumentor/timer.h"
#include "io/csv_writer_2d.h"
#include "pe/poisson/poisson_solver2d.h"

#include <string>

void ConcatPoissonSolver2D::set_parameter(int in_m, double in_tol, int in_maxIter)
{
    m       = in_m;
    tol     = in_tol;
    maxIter = in_maxIter;
}

ConcatPoissonSolver2D::ConcatPoissonSolver2D(Variable2D* in_variable)
    : variable(in_variable)
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

    EnvironmentConfig& env_cfg = EnvironmentConfig::Get();

    showGmresRes = env_cfg.showGmresRes;
}

void ConcatPoissonSolver2D::construct_solver_map_at_domain(Domain2DUniform* domain)
{
    if (tree_map[domain].size() > 0)
    {
        solver_map[domain] = new GMRESSolver2D(domain, m, tol, maxIter);
        auto* gmres        = static_cast<GMRESSolver2D*>(solver_map[domain]);
        gmres->set_solver(new PoissonSolver2D(domain, variable));
        gmres->schur_mat_construct(tree_map[domain], solver_map);
    }
    else
    {
        solver_map[domain] = new PoissonSolver2D(domain, variable);
    }
}

void ConcatPoissonSolver2D::construct_solver_map()
{
    SCOPE_TIMER("[Concat] Schur complement build total = ",
                TimeRecordType::None,
                EnvironmentConfig::Get().track_pe_construct_time);

    for (auto it = hierarchical_solve_levels.rbegin(); it != hierarchical_solve_levels.rend(); ++it)
        for (Domain2DUniform* domain : *it)
            construct_solver_map_at_domain(domain);

    construct_solver_map_at_domain(tree_root);
}

void ConcatPoissonSolver2D::solve()
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
            SCOPE_TIMER(env_cfg.pe_solve_total_name, TimeRecordType::Accumulate, false);
            solver_map[domain]->solve(*temp_fields[domain]);
        }
    }

    // Root equation
    for (auto& [location, child_domain] : tree_map[tree_root])
        field_map[tree_root]->bond_add(location, -1., *temp_fields[child_domain]);
    {
        SCOPE_TIMER(env_cfg.pe_solve_total_name, TimeRecordType::Accumulate, false);
        solver_map[tree_root]->solve(*field_map[tree_root]);
    }

    // Branch equations
    for (auto level : hierarchical_solve_levels)
    {
        for (Domain2DUniform* domain : level)
            field_map[domain]->bond_add(parent_map[domain].first, -1., *field_map[parent_map[domain].second]);
        for (Domain2DUniform* domain : level)
        {
            SCOPE_TIMER(env_cfg.pe_solve_total_name, TimeRecordType::Accumulate, false);
            solver_map[domain]->solve(*field_map[domain]);
        }
    }
}

void ConcatPoissonSolver2D::boundary_assembly()
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

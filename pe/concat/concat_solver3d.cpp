#include "concat_solver3d.h"
#include "gmres_solver3d.h"
#include "pe/poisson/poisson_solver3d.h"

void ConcatPoissonSolver3D::set_parameter(int in_m, double in_tol, int in_maxIter)
{
    const bool parameter_changed = (m != in_m) || (tol != in_tol) || (maxIter != in_maxIter);

    m       = in_m;
    tol     = in_tol;
    maxIter = in_maxIter;

    if (parameter_changed && !solver_map.empty())
    {
        clear_solver_map();
        construct_solver_map();
    }
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

    clear_solver_map();
}

void ConcatPoissonSolver3D::clear_solver_map()
{
    for (auto& [domain, solver] : solver_map)
        delete solver;

    solver_map.clear();
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
    SCOPE_TIMER("[Concat] Schur complement build total = ",
                TimeRecordType::None,
                EnvironmentConfig::Get().track_pe_construct_time);

    for (auto it = hierarchical_solve_levels.rbegin(); it != hierarchical_solve_levels.rend(); ++it)
        for (Domain3DUniform* domain : *it)
            construct_solver_map_at_domain(domain);

    construct_solver_map_at_domain(tree_root);
}

void ConcatPoissonSolver3D::solve()
{
    EnvironmentConfig& env_cfg = EnvironmentConfig::Get();

    SCOPE_TIMER("ConcatPoissonSolver3D::solve", TimeRecordType::None, env_cfg.track_pe_solve_detail_time);

    // Boundary
    boundary_assembly();

    // *hx*hx for each field
    for (auto kv : field_map)
    {
        auto  domain = kv.first;
        auto& f      = *kv.second;
        f *= domain->hx * domain->hx;
    }

    // XPositivehand construction (zneg-zpos pass)
    if (env_cfg.showCurrentStep)
        std::cout << "[Concat3D] Bottom-zpos pass: start" << std::endl;
    for (auto it = hierarchical_solve_levels.rbegin(); it != hierarchical_solve_levels.rend(); ++it)
    {
        auto level = *it;
        for (Domain3DUniform* domain : level)
        {
            if (env_cfg.showCurrentStep)
                std::cout << "[Concat3D] Bottom-zpos pass: solve domain " << domain->name << std::endl;
            (*temp_fields[domain]) = (*field_map[domain]);
            for (auto& [location, child_domain] : tree_map[domain])
            {
                temp_fields[domain]->bond_add(location, -1., *temp_fields[child_domain]);
                field_map[domain]->bond_add(location, -1., *temp_fields[child_domain]);
            }
        }
        for (Domain3DUniform* domain : level)
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

    for (auto level : hierarchical_solve_levels)
    {
        for (Domain3DUniform* domain : level)
            field_map[domain]->bond_add(parent_map[domain].first, -1., *field_map[parent_map[domain].second]);
        for (Domain3DUniform* domain : level)
        {
            SCOPE_TIMER(env_cfg.pe_solve_total_name, TimeRecordType::Accumulate, false);
            solver_map[domain]->solve(*field_map[domain]);
        }
    }
}

void ConcatPoissonSolver3D::boundary_assembly()
{
    // Apply boundary conditions to all domains in the geometry
    for (auto kv : field_map)
    {
        Domain3DUniform* domain = kv.first;
        field3&          f      = *kv.second;

        auto& var_value_map = variable->buffer_map[domain];
        auto& var_type_map  = variable->boundary_type_map[domain];

        int    nx = domain->get_nx();
        int    ny = domain->get_ny();
        int    nz = domain->get_nz();
        double hx = domain->get_hx();
        double hy = domain->get_hy();
        double hz = domain->get_hz();

        PDEBoundaryType boundary_type_xneg = var_type_map[LocationType::XNegative];
        PDEBoundaryType boundary_type_xpos = var_type_map[LocationType::XPositive];
        PDEBoundaryType boundary_type_yneg = var_type_map[LocationType::YNegative];
        PDEBoundaryType boundary_type_ypos = var_type_map[LocationType::YPositive];
        PDEBoundaryType boundary_type_zneg = var_type_map[LocationType::ZNegative];
        PDEBoundaryType boundary_type_zpos = var_type_map[LocationType::ZPositive];

        if (boundary_type_xneg == PDEBoundaryType::Dirichlet)
        {
            field2& boundary_value = *var_value_map[LocationType::XNegative];
            f.xneg_bond_add(-1.0 / hx / hx, boundary_value);
        }
        if (boundary_type_xpos == PDEBoundaryType::Dirichlet)
        {
            field2& boundary_value = *var_value_map[LocationType::XPositive];
            f.xpos_bond_add(-1.0 / hx / hx, boundary_value);
        }

        if (boundary_type_yneg == PDEBoundaryType::Dirichlet)
        {
            field2& boundary_value = *var_value_map[LocationType::YNegative];
            f.yneg_bond_add(-1.0 / hy / hy, boundary_value);
        }
        if (boundary_type_ypos == PDEBoundaryType::Dirichlet)
        {
            field2& boundary_value = *var_value_map[LocationType::YPositive];
            f.back_bond_add(-1.0 / hy / hy, boundary_value);
        }

        if (boundary_type_zneg == PDEBoundaryType::Dirichlet)
        {
            field2& boundary_value = *var_value_map[LocationType::ZNegative];
            f.zneg_bond_add(-1.0 / hz / hz, boundary_value);
        }
        if (boundary_type_zpos == PDEBoundaryType::Dirichlet)
        {
            field2& boundary_value = *var_value_map[LocationType::ZPositive];
            f.up_bond_add(-1.0 / hz / hz, boundary_value);
        }

        if (boundary_type_xneg == PDEBoundaryType::Neumann)
        {
            field2& boundary_value = *var_value_map[LocationType::XNegative];
            f.xneg_bond_add(1.0 / hx, boundary_value);
        }
        if (boundary_type_xpos == PDEBoundaryType::Neumann)
        {
            field2& boundary_value = *var_value_map[LocationType::XPositive];
            f.xpos_bond_add(-1.0 / hx, boundary_value);
        }
        if (boundary_type_yneg == PDEBoundaryType::Neumann)
        {
            field2& boundary_value = *var_value_map[LocationType::YNegative];
            f.yneg_bond_add(1.0 / hy, boundary_value);
        }
        if (boundary_type_ypos == PDEBoundaryType::Neumann)
        {
            field2& boundary_value = *var_value_map[LocationType::YPositive];
            f.back_bond_add(-1.0 / hy, boundary_value);
        }
        if (boundary_type_zneg == PDEBoundaryType::Neumann)
        {
            field2& boundary_value = *var_value_map[LocationType::ZNegative];
            f.zneg_bond_add(1.0 / hz, boundary_value);
        }
        if (boundary_type_zpos == PDEBoundaryType::Neumann)
        {
            field2& boundary_value = *var_value_map[LocationType::ZPositive];
            f.up_bond_add(-1.0 / hz, boundary_value);
        }
    }
}

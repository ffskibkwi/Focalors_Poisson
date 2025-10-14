#include "concat_solver2d.h"

ConcatSolver2D::ConcatSolver2D(Variable &in_variable, EnvironmentConfig* in_env_config)
    : variable(in_variable)
    , env_config(in_env_config)
{
    //config load
    if (in_env_config)
    {
        showGmresRes = in_env_config->showGmresRes;
    }

    // geometry double check
    if (variable.geometry == nullptr)
        throw std::runtime_error("ConcatSolver2D: variable.geometry is null");
    if (!variable.geometry->is_checked)
        variable.geometry->check();
    if (variable.geometry->tree_root == nullptr || variable.geometry->tree_map.empty())
        variable.geometry->solve_prepare();

    tree_root = variable.geometry->tree_root;
    tree_map  = variable.geometry->tree_map;
    parent_map= variable.geometry->parent_map;
    field_map = variable.field_map;

    specify_solve_order();
    construct_solver_map();
    
    //Construct the temp field for each domain
    for (auto &[domain, field] : field_map)
    {
        if (domain != tree_root)
            temp_fields[domain] = new field2(field->get_nx(), field->get_ny(), field->get_name() + "_temp");
    }
}

ConcatSolver2D::~ConcatSolver2D()
{
    for (auto &[domain, temp_field] : temp_fields)
        delete temp_field;
    for (auto &domain : solve_order)
        delete solver_map[domain];
}

void ConcatSolver2D::specify_solve_order()
{
    // Solve order arrangement
    std::queue<Domain2DUniform*> q;
    q.push(tree_root);
    while(!q.empty())
    {
        Domain2DUniform* current = q.front();
        q.pop();
        if (current != tree_root)
            solve_order.insert(solve_order.begin(), current);
        if (tree_map.count(current))
        {
            for (auto &kv : tree_map[current])
                q.push(kv.second);
        }
    }
}

void ConcatSolver2D::construct_solver_map()
{
    //Construct solvers (for non-root domains)
    for (auto &domain : solve_order)
    {
        if (tree_map[domain].size() > 0)
        {
            solver_map[domain] = new GMRESSolver2D(domain, m, tol, maxIter, env_config);
            static_cast<GMRESSolver2D*>(solver_map[domain])->schur_mat_construct(tree_map[domain], solver_map); //Here use RTTI
        }    
        else
        {
            solver_map[domain] = new PoissonSolver2D(domain, env_config);
        } 
    }

    //Construct solver for root domain
    if (tree_map[tree_root].size() > 0)
    {
        solver_map[tree_root] = new GMRESSolver2D(tree_root, m, tol, maxIter, env_config);
        static_cast<GMRESSolver2D*>(solver_map[tree_root])->schur_mat_construct(tree_map[tree_root], solver_map); //Here use RTTI
    }    
    else
    {
        solver_map[tree_root] = new PoissonSolver2D(tree_root, env_config);
    } 
}

void ConcatSolver2D::solve()
{
    //Righthand construction
    for (auto &domain : solve_order)
    {
        (*temp_fields[domain]) = (*field_map[domain]);
        for (auto &[location, child_domain] : tree_map[domain])
        {
            temp_fields[domain]->bond_add(location, -1., *temp_fields[child_domain]);
            field_map[domain]->bond_add(location, -1., *temp_fields[child_domain]);
        }
        solver_map[domain]->solve(*temp_fields[domain]);
    }
    
    //Root equation
    for (auto &[location, child_domain] : tree_map[tree_root])
        field_map[tree_root]->bond_add(location, -1., *temp_fields[child_domain]);
    solver_map[tree_root]->solve(*field_map[tree_root]);

    //Branch equations
    for (auto it = solve_order.rbegin(); it != solve_order.rend(); ++it)
    {
        Domain2DUniform* d = *it;
        field_map[d]->bond_add(parent_map[d].first, -1., *field_map[parent_map[d].second]);
        solver_map[d]->solve(*field_map[d]);
    }
}
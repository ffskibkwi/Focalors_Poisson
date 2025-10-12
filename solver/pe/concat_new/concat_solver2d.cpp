#include "concat_solver2d.h"

ConcatSolver2D::ConcatSolver2D(Variable &in_variable)
    : variable(in_variable)
{
    specify_solve_order();
    construct_solver_map();
    
    //Construct the temp field for each domain
    for (auto &[domain, field] : variable.field_map)
        temp_fields[&domain] = new field2(field.get_nx(), field.get_ny(), field.get_name() + "_temp");
}

~ConcatSolver2D()
{
    for (auto &[domain, temp_field] : temp_fields)
        delete temp_field;
    for (auto &domain : solve_order)
        delete solver_map[domain];
}

void ConcatSolver2D::specify_solve_order()
{
    GeometryTreeNode<Domain2DUniform>* root_domain = variable.geometry.tree_root;
    
    // Solve order arrangement
    std::queue<Domain2DUniform*> q;
    q.push(root_domain->domain);
    while(!q.empty())
    {
        Domain2DUniform* current = q.front();
        q.pop();
        solve_order.insert(solve_order.begin(), current);
        for (auto &child : current->children)
        {
            q.push(child->domain);  //In this version, assume there is no ring, so do not check is child already in q
        }
    }
}

void ConcatSolver2D::construct_solver_map()
{
    //Construct solvers
    for (auto &domain : solve_order)
    {
        if (variable.geometry.adjacency[domain].size() > 0)
        {
            solver_map[domain] = new GMRESSolver2D(domain, m, tol, maxIter);
            solver_map[domain]->schur_mat_construct(variable.geometry.adjacency[domain], solver_map);
        }    
        else
        {
            solver_map[domain] = new PoissonSolver2D(domain);
        } 
    }
}

void ConcatSolver2D::solve()
{
    for (auto &domain : solve_order)
    {
        (*temp_fields[domain]) = variable.field_map[domain];
        for (auto &[location, child_domain] : variable.geometry.tree_map[domain])
        {
            temp_fields[domain]->bond_add(location, -1., *temp_fields[child_domain]);
        }
        solver_map[domain]->solve(*temp_fields[domain]);
    }
}
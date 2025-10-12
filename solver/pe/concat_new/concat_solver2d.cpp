#include "concat_solver2d.h"

ConcatSolver2D::ConcatSolver2D(Variable &in_variable)
    : variable(in_variable)
{
    init();
}



void ConcatSolver2D::construct_solve_stack()
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


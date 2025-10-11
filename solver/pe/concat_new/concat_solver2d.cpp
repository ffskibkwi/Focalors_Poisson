#include "concat_solver2d.h"

ConcatSolver2D::ConcatSolver2D(Variable &in_variable)
    : variable(in_variable)
{
    init();
}



void ConcatSolver2D::construct_solve_stack()
{
    GeometryTreeNode* root_domain = variable.geometry.tree_root;

    // for (auto &[domain, field] : variable.field_map)
    // {
    //     solve_queue.push(&domain);
    // }
}


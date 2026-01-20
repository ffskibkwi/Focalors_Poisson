#include "concat_solver3d.h"

void ConcatPoissonSolver3D::set_parameter(int in_m, double in_tol, int in_maxIter)
{
    m = in_m;
    tol = in_tol;
    maxIter = in_maxIter;
}

ConcatPoissonSolver3D::ConcatPoissonSolver3D(Variable3D* in_variable, EnvironmentConfig* in_env_config)
    : variable(in_variable)
    , env_config(in_env_config)
{
    //config load
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

    tree_root = variable->geometry->tree_root;
    tree_map  = variable->geometry->tree_map;
    parent_map= variable->geometry->parent_map;
    field_map = variable->field_map;

    specify_solve_order();
    construct_solver_map();
    
    //Construct the temp field for each domain
    for (auto &[domain, field] : field_map)
    {
        if (domain != tree_root)
            temp_fields[domain] = new field3(field->get_nx(), field->get_ny(), field->get_nz(), field->get_name() + "_temp");
    }
}

ConcatPoissonSolver3D::~ConcatPoissonSolver3D()
{
    for (auto &[domain, temp_field] : temp_fields)
        delete temp_field;
    for (auto &domain : solve_order)
        delete solver_map[domain];
    if (tree_root)
        delete solver_map[tree_root];
}

void ConcatPoissonSolver3D::specify_solve_order()
{
    // Solve order arrangement
    std::queue<Domain3DUniform*> q;
    q.push(tree_root);
    while(!q.empty())
    {
        Domain3DUniform* current = q.front();
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

void ConcatPoissonSolver3D::construct_solver_map()
{
    //Construct solvers (for non-root domains)
    for (auto &domain : solve_order)
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

    //Construct solver for root domain
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

void ConcatPoissonSolver3D::bond_add_3d(field3& target, LocationType location, double k, const field3& source)
{
    int target_nx = target.get_nx();
    int target_ny = target.get_ny();
    int target_nz = target.get_nz();
    int source_nx = source.get_nx();
    int source_ny = source.get_ny();
    int source_nz = source.get_nz();

    switch (location)
    {
        case LocationType::Left:
            // Left face: x=0 plane
            for (int j = 0; j < target_ny; j++)
                for (int k = 0; k < target_nz; k++)
                    target(0, j, k) += k * source(source_nx - 1, j, k);
            break;
        case LocationType::Right:
            // Right face: x=nx-1 plane
            for (int j = 0; j < target_ny; j++)
                for (int k = 0; k < target_nz; k++)
                    target(target_nx - 1, j, k) += k * source(0, j, k);
            break;
        case LocationType::Front:
            // Front face: y=0 plane
            for (int i = 0; i < target_nx; i++)
                for (int k = 0; k < target_nz; k++)
                    target(i, 0, k) += k * source(i, source_ny - 1, k);
            break;
        case LocationType::Back:
            // Back face: y=ny-1 plane
            for (int i = 0; i < target_nx; i++)
                for (int k = 0; k < target_nz; k++)
                    target(i, target_ny - 1, k) += k * source(i, 0, k);
            break;
        case LocationType::Down:
            // Down face: z=0 plane
            for (int i = 0; i < target_nx; i++)
                for (int j = 0; j < target_ny; j++)
                    target(i, j, 0) += k * source(i, j, source_nz - 1);
            break;
        case LocationType::Up:
            // Up face: z=nz-1 plane
            for (int i = 0; i < target_nx; i++)
                for (int j = 0; j < target_ny; j++)
                    target(i, j, target_nz - 1) += k * source(i, j, 0);
            break;
        default:
            throw std::invalid_argument("Invalid location type for 3D bond_add");
    }
}

void ConcatPoissonSolver3D::solve()
{
    //Righthand construction (bottom-up pass)
    for (auto &domain : solve_order)
    {
        (*temp_fields[domain]) = (*field_map[domain]);
        for (auto &[location, child_domain] : tree_map[domain])
        {
            bond_add_3d(*temp_fields[domain], location, -1., *temp_fields[child_domain]);
            bond_add_3d(*field_map[domain], location, -1., *temp_fields[child_domain]);
        }
        if (env_config && env_config->showCurrentStep)
        {
            double s_pre = temp_fields[domain]->sum();
            std::cout << "[Concat3D] Domain " << domain->name << " temp sum before solve=" << s_pre << std::endl;
        }
        solver_map[domain]->solve(*temp_fields[domain]);
        if (env_config && env_config->showCurrentStep)
        {
            double s_post = temp_fields[domain]->sum();
            std::cout << "[Concat3D] Domain " << domain->name << " temp sum after solve=" << s_post << std::endl;
        }
    }
    
    //Root equation
    for (auto &[location, child_domain] : tree_map[tree_root])
        bond_add_3d(*field_map[tree_root], location, -1., *temp_fields[child_domain]);
    if (env_config && env_config->showCurrentStep)
    {
        double s_root_pre = field_map[tree_root]->sum();
        std::cout << "[Concat3D] Root domain " << tree_root->name << " field sum before solve=" << s_root_pre << std::endl;
    }
    solver_map[tree_root]->solve(*field_map[tree_root]);
    if (env_config && env_config->showCurrentStep)
    {
        double s_root_post = field_map[tree_root]->sum();
        std::cout << "[Concat3D] Root domain " << tree_root->name << " field sum after solve=" << s_root_post << std::endl;
    }

    //Branch equations (top-down pass)
    for (auto it = solve_order.rbegin(); it != solve_order.rend(); ++it)
    {
        Domain3DUniform* d = *it;
        bond_add_3d(*field_map[d], parent_map[d].first, -1., *field_map[parent_map[d].second]);
        if (env_config && env_config->showCurrentStep)
        {
            double s_pre = field_map[d]->sum();
            std::cout << "[Concat3D] Branch domain " << d->name << " field sum before solve=" << s_pre << std::endl;
        }
        solver_map[d]->solve(*field_map[d]);
        if (env_config && env_config->showCurrentStep)
        {
            double s_post = field_map[d]->sum();
            std::cout << "[Concat3D] Branch domain " << d->name << " field sum after solve=" << s_post << std::endl;
        }
    }
}
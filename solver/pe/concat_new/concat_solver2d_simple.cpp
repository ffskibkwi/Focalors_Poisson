#include "concat_solver2d.h"
#include "pe/poisson_new/poisson_solver2d.h"

ConcatSolver2D_Simple::ConcatSolver2D_Simple(Variable &in_variable)
    : variable(in_variable)
{
    init();
}


void ConcatSolver2D_Simple::init()
{
    geo_main_domain = variable.geometry.main_domain;
    geo_main_neighbour_domains = variable.geometry->adjacency[geo_main_domain];
    
    for (auto &[domain, field] : variable.field_map)
    {
        //Construct the Poisson solver for each domain
        pe_solver_map[&domain] = new PoissonSolver2D(domain);
        //Construct the temp field for each domain
        temp_fields[&domain] = new field2(field.get_nx(), field.get_ny(), field.get_name() + "_temp");
    }

    //In current version, the geometry with only one main domain is supported
    //Thus the Schur matrix construction has a "core", i.e., the main domain
    for (auto &[LocationType, domain] : geo_main_neighbour_domains)
    {
        //Construct the Schur matrix for each neibour domain of main domain
        switch (LocationType)
        {
            case LocationType::Left:
            {
                schur_mat_map[&domain] = new Schur_mat_left(geo_main_domain, domain);
                schur_mat_map[&domain]->construct(pe_solver_map[&domain]);
                schur_mat_vec.push_back(schur_mat_map[&domain]);
            }  
            break;
            case LocationType::Right:
            {
                schur_mat_map[&domain] = new Schur_mat_right(geo_main_domain, domain);
                schur_mat_map[&domain]->construct(pe_solver_map[&domain]);
                schur_mat_vec.push_back(schur_mat_map[&domain]);
            }
            break;
            case LocationType::Up:
            {
                schur_mat_map[&domain] = new Schur_mat_up(geo_main_domain, domain);
                schur_mat_map[&domain]->construct(pe_solver_map[&domain]);
                schur_mat_vec.push_back(schur_mat_map[&domain]);
            }
            break;
            case LocationType::Down:
            {
                schur_mat_map[&domain] = new Schur_mat_down(geo_main_domain, domain);
                schur_mat_map[&domain]->construct(pe_solver_map[&domain]);
                schur_mat_vec.push_back(schur_mat_map[&domain]);
            }
            break;
            default:
                throw std::invalid_argument("Invalid location type");
        }
    }
}


void ConcatSolver2D_Simple::solve()
{
    // for (auto &[domain, field] : variable.field_map)     // It is error, should not including main field
    for (auto &[LocationType, domain] : geo_main_neighbour_domains)
    {
        (*temp_fields[&domain]) = (*domain);
        pe_solver_map[&domain]->solve(*temp_fields[&domain]);
    }
    
    for (auto &[LocationType, domain] : geo_main_neighbour_domains)
    {
        //Construct the Schur matrix for each neibour domain of main domain
        switch (LocationType)
        {
            case LocationType::Left:
                geo_main_domain->left_bond_add(-1., *temp_fields[&domain]);
            break;
            case LocationType::Right:
                geo_main_domain->right_bond_add(-1., *temp_fields[&domain]);
            break;
            case LocationType::Up:
                geo_main_domain->up_bond_add(-1., *temp_fields[&domain]);
            break;
            case LocationType::Down:
                geo_main_domain->down_bond_add(-1., *temp_fields[&domain]);
            break;
            default:
                throw std::invalid_argument("Invalid location type");
        }
    }

    (*geo_main_domain) = gmres(*geo_main_domain, *geo_main_domain, schur_mat_vec, pe_solver_map[geo_main_domain], gmres_m, gmres_tol, 1000, resVec);

    for (auto &[LocationType, domain] : geo_main_neighbour_domains)
    {
        //Construct the Schur matrix for each neibour domain of main domain
        switch (LocationType)
        {
            case LocationType::Left:
                domain->right_bond_add(-1., *geo_main_domain);
            break;
            case LocationType::Right:
                domain->left_bond_add(-1., *geo_main_domain);
            break;
            case LocationType::Up:
                domain->down_bond_add(-1., *geo_main_domain);
            break;
            case LocationType::Down:
                domain->up_bond_add(-1., *geo_main_domain);
            break;
            default:
                throw std::invalid_argument("Invalid location type");
        }
    }

    for (auto &[LocationType, domain] : geo_main_neighbour_domains)
        pe_solver_map[&domain]->solve(*temp_fields[&domain]);
}



/*

std::unordered_map<Domain2DUniform*, std::unordered_map<LocationType, Domain2DUniform*>> adjacency;

field2 p_1(nx_1, ny_1, "p_1");
field2 p_2(nx_2, ny_2, "p_2");
field2 p_3(nx_3, ny_3, "p_3");
field2 p_4(nx_4, ny_4, "p_4");
field2 p_5(nx_5, ny_5, "p_5");

field2 p_temp_1(nx_1, ny_1, "p_temp_1");
field2 p_temp_3(nx_3, ny_3, "p_temp_3");
field2 p_temp_4(nx_4, ny_4, "p_temp_4");
field2 p_temp_5(nx_5, ny_5, "p_temp_5");

PoissonSolver2D<PDEBoundaryType::Dirichlet,
                PDEBoundaryType::Dirichlet,
                PDEBoundaryType::Neumann,
                PDEBoundaryType::Neumann>
    pe_solver_1(p_1.get_nx(), p_1.get_ny());
PoissonSolver2D<PDEBoundaryType::Dirichlet,
                PDEBoundaryType::Dirichlet,
                PDEBoundaryType::Dirichlet,
                PDEBoundaryType::Dirichlet>
    pe_solver_2(p_2.get_nx(), p_2.get_ny());
PoissonSolver2D<PDEBoundaryType::Dirichlet,
                PDEBoundaryType::Dirichlet,
                PDEBoundaryType::Neumann,
                PDEBoundaryType::Neumann>
    pe_solver_3(p_3.get_nx(), p_3.get_ny());
PoissonSolver2D<PDEBoundaryType::Neumann,
                PDEBoundaryType::Neumann,
                PDEBoundaryType::Dirichlet,
                PDEBoundaryType::Dirichlet>
    pe_solver_4(p_4.get_nx(), p_4.get_ny());
PoissonSolver2D<PDEBoundaryType::Neumann,
                PDEBoundaryType::Neumann,
                PDEBoundaryType::Dirichlet,
                PDEBoundaryType::Dirichlet>
    pe_solver_5(p_5.get_nx(), p_5.get_ny());

Shur_mat_left  S_21(p_2, p_1);
Shur_mat_right S_23(p_2, p_3);
Shur_mat_down  S_24(p_2, p_4);
Shur_mat_up    S_25(p_2, p_5);

S_21.construct(pe_solver_1);
S_23.construct(pe_solver_3);
S_24.construct(pe_solver_4);
S_25.construct(pe_solver_5);



//Solve
    p_temp_1 = p_1;
    p_temp_3 = p_3;
    p_temp_4 = p_4;
    p_temp_5 = p_5;

    pe_solver_1.solve(p_temp_1);
    pe_solver_3.solve(p_temp_3);
    pe_solver_4.solve(p_temp_4);
    pe_solver_5.solve(p_temp_5);

    p_2.left_bond_add(-1., p_temp_1);
    p_2.right_bond_add(-1., p_temp_3);
    p_2.down_bond_add(-1., p_temp_4);
    p_2.up_bond_add(-1., p_temp_5);
    pe_solver_2.solve(p_2);

    std::vector<double> resVec;
    // std::vector<Shur_mat*> S = {&S_21, &S_23, &S_24, &S_25};
    // p_2 = gmres(p_2, p_2, S, solver_2, 5, 1.e-7, 100, resVec);
    // p_2 = gmres(p_2, p_2, S_21, S_23, solver_2, 3, 1.e-7, 100, resVec);

    // std::cout << "acc_rank:" << acc_rank << ", gmres_m:" << gmres_m << ", result:";
    std::cout << acc_rank << ", " << gmres_m << ", "; //for csv output

    p_2 = gmres(p_2, p_2, S_21, S_23, S_24, S_25, pe_solver_2, gmres_m, gmres_tol, 1000, resVec);
    for (size_t i = 0; i < resVec.size(); ++i)
    {
        std::cout << resVec[i];
        if (i < resVec.size() - 1)
        {
            std::cout << ", ";
        }
    }
    std::cout << std::endl;
    p_1.right_bond_add(-1., p_2);
    p_3.left_bond_add(-1., p_2);
    p_4.up_bond_add(-1., p_2);
    p_5.down_bond_add(-1., p_2);
    pe_solver_1.solve(p_1);
    pe_solver_3.solve(p_3);
    pe_solver_4.solve(p_4);
    pe_solver_5.solve(p_5);

*/
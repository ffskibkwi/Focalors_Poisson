#include "core/boundary/boundary_type.h"
#include "io/csv_writer_2d.h"
#include "pe/concat/Shur_mat.hpp"
#include "pe/concat/Shur_solver.hpp"
#include "pe/concat/gmres.h"
#include "pe/poisson/poisson_solver_2d.hpp"

#include <functional>

int main(int argc, char* argv[])
{
    double acc_rank  = 1.0;
    int    base_nx_1 = 4;
    int    base_nx_2 = 8;
    int    base_nx_3 = 16;
    int    base_ny_4 = 8;
    int    base_ny_2 = 16;
    int    base_ny_5 = 4;

    int nx_1 = static_cast<int>(std::round(base_nx_1 * acc_rank));
    int nx_2 = static_cast<int>(std::round(base_nx_2 * acc_rank));
    int nx_3 = static_cast<int>(std::round(base_nx_3 * acc_rank));
    int ny_4 = static_cast<int>(std::round(base_ny_4 * acc_rank));
    int ny_2 = static_cast<int>(std::round(base_ny_2 * acc_rank));
    int ny_5 = static_cast<int>(std::round(base_ny_5 * acc_rank));
    int ny_1 = ny_2;
    int ny_3 = ny_2;
    int nx_4 = nx_2;
    int nx_5 = nx_2;

    double H  = 1. / (nx_1 + nx_2 + nx_3 + 1.);
    double L1 = (nx_1 + .5) * H;
    double L2 = nx_2 * H;
    double L3 = (nx_3 + .5) * H;
    double H4 = (ny_4 + .5) * H;
    double H2 = ny_2 * H;
    double H5 = (ny_5 + .5) * H;

    int    gmres_m   = 3;
    double gmres_tol = 1.e-10;

    // FluidContext2DNonExtended fluid_context_1;
    // FluidContext2DNonExtended fluid_context_2;
    // FluidContext2DNonExtended fluid_context_3;
    // FluidContext2DNonExtended fluid_context_4;
    // FluidContext2DNonExtended fluid_context_5;
    // fluid_context_1.init(nx_1, ny_1, H);
    // fluid_context_2.init(nx_2, ny_2, H);
    // fluid_context_3.init(nx_3, ny_3, H);
    // fluid_context_4.init(nx_4, ny_4, H);
    // fluid_context_5.init(nx_5, ny_5, H);

    // auto& p_1 = fluid_context_1.add_field(nx_1, ny_1, "p");
    // auto& p_2 = fluid_context_2.add_field(nx_2, ny_2, "p");
    // auto& p_3 = fluid_context_3.add_field(nx_3, ny_3, "p");
    // auto& p_4 = fluid_context_4.add_field(nx_4, ny_4, "p");
    // auto& p_5 = fluid_context_5.add_field(nx_5, ny_5, "p");

    // // temporary perssure for pre calc
    // auto& p_temp_1 = fluid_context_1.add_field(nx_1, ny_1, "p_temp");
    // auto& p_temp_3 = fluid_context_3.add_field(nx_3, ny_3, "p_temp");
    // auto& p_temp_4 = fluid_context_4.add_field(nx_4, ny_4, "p_temp");
    // auto& p_temp_5 = fluid_context_5.add_field(nx_5, ny_5, "p_temp");

    field2 p_1(nx_1, ny_1, "p_1");
    field2 p_2(nx_2, ny_2, "p_2");
    field2 p_3(nx_3, ny_3, "p_3");
    field2 p_4(nx_4, ny_4, "p_4");
    field2 p_5(nx_5, ny_5, "p_5");

    field2 p_temp_1(nx_1, ny_1, "p_temp_1");
    field2 p_temp_3(nx_3, ny_3, "p_temp_3");
    field2 p_temp_4(nx_4, ny_4, "p_temp_4");
    field2 p_temp_5(nx_5, ny_5, "p_temp_5");

    PoissonSolver2D<FluidBoundaryType::Dirichlet,
                    FluidBoundaryType::Dirichlet,
                    FluidBoundaryType::Neumann,
                    FluidBoundaryType::Neumann>
        pe_solver_1(p_1.get_nx(), p_1.get_ny());
    PoissonSolver2D<FluidBoundaryType::Dirichlet,
                    FluidBoundaryType::Dirichlet,
                    FluidBoundaryType::Dirichlet,
                    FluidBoundaryType::Dirichlet>
        pe_solver_2(p_2.get_nx(), p_2.get_ny());
    PoissonSolver2D<FluidBoundaryType::Dirichlet,
                    FluidBoundaryType::Dirichlet,
                    FluidBoundaryType::Neumann,
                    FluidBoundaryType::Neumann>
        pe_solver_3(p_3.get_nx(), p_3.get_ny());
    PoissonSolver2D<FluidBoundaryType::Neumann,
                    FluidBoundaryType::Neumann,
                    FluidBoundaryType::Dirichlet,
                    FluidBoundaryType::Dirichlet>
        pe_solver_4(p_4.get_nx(), p_4.get_ny());
    PoissonSolver2D<FluidBoundaryType::Neumann,
                    FluidBoundaryType::Neumann,
                    FluidBoundaryType::Dirichlet,
                    FluidBoundaryType::Dirichlet>
        pe_solver_5(p_5.get_nx(), p_5.get_ny());

    Shur_mat_left  S_21(p_2, p_1);
    Shur_mat_right S_23(p_2, p_3);
    Shur_mat_down  S_24(p_2, p_4);
    Shur_mat_up    S_25(p_2, p_5);

    S_21.construct(pe_solver_1);
    S_23.construct(pe_solver_3);
    S_24.construct(pe_solver_4);
    S_25.construct(pe_solver_5);

    auto set_value = [&](double (*op)(double, double, double)) {
        for (int i = 0; i < p_1.get_nx(); i++)
        {
            for (int j = 0; j < p_1.get_ny(); j++)
            {
                p_1(i, j) = op((i + 0.5) * H, (j + 0.5) * H, p_1(i, j));
            }
        }

        for (int i = 0; i < p_2.get_nx(); i++)
        {
            for (int j = 0; j < p_2.get_ny(); j++)
            {
                p_2(i, j) = op((i + nx_1 + 0.5) * H, (j + 0.5) * H, p_2(i, j));
            }
        }

        for (int i = 0; i < p_3.get_nx(); i++)
        {
            for (int j = 0; j < p_3.get_ny(); j++)
            {
                p_3(i, j) = op((i + nx_1 + 0.5) * H, (j - ny_3 + 0.5) * H, p_3(i, j));
            }
        }
    };

    // auto set_value_global = [&](double (*op)(double, double)) {
    auto set_value_global = [&](const std::function<double(double, double)>& op) {
        for (int i = 0; i < p_1.get_nx(); i++)
        {
            for (int j = 0; j < p_1.get_ny(); j++)
            {
                p_1(i, j) = op((i + 1.) * H, (j + ny_4 + 1.) * H);
            }
        }

        for (int i = 0; i < p_2.get_nx(); i++)
        {
            for (int j = 0; j < p_2.get_ny(); j++)
            {
                p_2(i, j) = op((i + nx_1 + 1.) * H, (j + ny_4 + 1.) * H);
            }
        }

        for (int i = 0; i < p_3.get_nx(); i++)
        {
            for (int j = 0; j < p_3.get_ny(); j++)
            {
                p_3(i, j) = op((i + nx_1 + nx_2 + 1.) * H, (j + ny_4 + 1.) * H);
            }
        }

        for (int i = 0; i < p_4.get_nx(); i++)
        {
            for (int j = 0; j < p_4.get_ny(); j++)
            {
                p_4(i, j) = op((i + nx_1 + 1.) * H, (j + 1.) * H);
            }
        }

        for (int i = 0; i < p_5.get_nx(); i++)
        {
            for (int j = 0; j < p_5.get_ny(); j++)
            {
                p_5(i, j) = op((i + nx_1 + 1.) * H, (j + ny_4 + ny_2 + 1.) * H);
            }
        }
    };

    // auto circle1 = [](double x, double y, double raw) {
    //     double cx = 20.0;
    //     double cy = 0.0;
    //     double r  = 0.1;
    //     if ((x - cx) * (x - cx) + (y - cy) * (y - cy) <= r * r)
    //     {
    //         return 1.0;
    //     }
    //     return raw;
    // };

    // set_value(circle1);

    // 示例解析源项 f_analy（未直接使用）；此处保持为注释以避免捕获错误
    // auto f_analy = [L1, L2, L3, H4, H2, H5, H](double x, double y) {
    //     double A1_x           = M_PI / (2.0 * L1);
    //     double A2_x           = (M_PI / L2 - A1_x) / (L1 + L2);
    //     double numerator_A3_x = -((M_PI / L2 + M_PI / (2.0 * L3))) / (L2 + L3) - A2_x;
    //     double A3_x           = numerator_A3_x / (L1 + L2 + L3);
    //     double xi2_x    = A2_x + A3_x * (x - L1 - L2);
    //     double xi1_x    = A1_x + xi2_x * (x - L1);
    //     double psix_val = x * xi1_x;
    //     double xi1d_x     = xi2_x + A3_x * (x - L1);
    //     double dpsix_val  = xi1_x + x * xi1d_x;
    //     double ddpsix_val = 2.0 * xi1d_x + 2.0 * A3_x * x;
    //     double A1_y           = M_PI / (2.0 * H4);
    //     double A2_y           = (M_PI / H2 - A1_y) / (H4 + H2);
    //     double numerator_A3_y = -((M_PI / H2 + M_PI / (2.0 * H5))) / (H2 + H5) - A2_y;
    //     double A3_y           = numerator_A3_y / (H4 + H2 + H5);
    //     double xi2_y    = A2_y + A3_y * (y - H4 - H2);
    //     double xi1_y    = A1_y + xi2_y * (y - H4);
    //     double psiy_val = y * xi1_y;
    //     double xi1d_y     = xi2_y + A3_y * (y - H4);
    //     double dpsiy_val  = xi1_y + y * xi1d_y;
    //     double ddpsiy_val = 2.0 * xi1d_y + 2.0 * A3_y * y;
    //     double p_ana_val   = std::sin(psix_val) * std::sin(psiy_val);
    //     double pdd_ana_val = ddpsix_val * std::cos(psix_val) * std::sin(psiy_val) +
    //                          ddpsiy_val * std::cos(psiy_val) * std::sin(psix_val) -
    //                          p_ana_val * (dpsix_val * dpsiy_val + dpsiy_val * dpsiy_val);
    //     return pdd_ana_val * H * H;
    // };

    // 使用常数函数初始化
    auto f_const = [](double x, double y) { return 0.0; };
    set_value_global(f_const);

    // Calculation
    // SCOPE_TIMER 宏未定义，注释掉
    // SCOPE_TIMER("PE Solver 2 Solve");

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

    // IO::field_to_csv(p_1, root_dir + "/p_1");
    // IO::field_to_csv(p_2, root_dir + "/p_2");
    // IO::field_to_csv(p_3, root_dir + "/p_3");
    // IO::field_to_csv(p_4, root_dir + "/p_4");
    // IO::field_to_csv(p_5, root_dir + "/p_5");

    return 0;
}

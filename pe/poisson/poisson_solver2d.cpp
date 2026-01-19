#include "poisson_solver2d.h"
#include "io/csv_writer_2d.h"
#include <string>

PoissonSolver2D::PoissonSolver2D(int                in_nx,
                                 int                in_ny,
                                 double             in_hx,
                                 double             in_hy,
                                 PDEBoundaryType    in_boundary_type_left,
                                 PDEBoundaryType    in_boundary_type_right,
                                 PDEBoundaryType    in_boundary_type_down,
                                 PDEBoundaryType    in_boundary_type_up,
                                 EnvironmentConfig* in_env_config)
    : nx(in_nx)
    , ny(in_ny)
    , hx(in_hx)
    , hy(in_hy)
    , boundary_type_left(in_boundary_type_left)
    , boundary_type_right(in_boundary_type_right)
    , boundary_type_down(in_boundary_type_down)
    , boundary_type_up(in_boundary_type_up)
{
    env_config = in_env_config;
    init();
}

PoissonSolver2D::PoissonSolver2D(Domain2DUniform* in_domain, Variable* in_variable, EnvironmentConfig* in_env_config)
    : domain(in_domain)
    , var(in_variable)
    , nx(in_domain->nx)
    , ny(in_domain->ny)
    , hx(in_domain->hx)
    , hy(in_domain->hy)
{
    env_config = in_env_config;

    // For PoissonSolver, var is only use to define the boundary type, the boundary value is not defined here but in
    // ConcatSolver
    boundary_type_left  = var->boundary_type_map[domain][LocationType::Left];
    boundary_type_right = var->boundary_type_map[domain][LocationType::Right];
    boundary_type_down  = var->boundary_type_map[domain][LocationType::Down];
    boundary_type_up    = var->boundary_type_map[domain][LocationType::Up];

    init();
}

void PoissonSolver2D::init()
{
    buffer.init(
        nx, ny, "buffer"); // Here init a new field for buffer, but it can be replaced by another existed buffer field

    auto isDirLike = [](PDEBoundaryType t) {
        return t == PDEBoundaryType::Dirichlet || t == PDEBoundaryType::Adjacented;
    };

    if (boundary_type_down == PDEBoundaryType::Periodic && boundary_type_up == PDEBoundaryType::Periodic)
    {
        poisson_fft_y = new PoissonFFT2D_PP();
    }
    else if (boundary_type_down == PDEBoundaryType::Neumann && boundary_type_up == PDEBoundaryType::Neumann)
    {
        poisson_fft_y = new PoissonFFT2D_NN();
    }
    else if (isDirLike(boundary_type_down) && isDirLike(boundary_type_up))
    {
        poisson_fft_y = new PoissonFFT2D_DD();
    }
    else if (isDirLike(boundary_type_down) && boundary_type_up == PDEBoundaryType::Neumann)
    {
        poisson_fft_y = new PoissonFFT2D_DN();
    }
    else if (boundary_type_down == PDEBoundaryType::Neumann && isDirLike(boundary_type_up))
    {
        poisson_fft_y = new PoissonFFT2D_ND();
    }
    poisson_fft_y->init(nx, ny);

    lambda_y = new double[ny];
    x_diag   = new double[ny];

    cal_lambda();

    for (int j = 0; j < ny; j++)
    {
        x_diag[j] = -2.0 + hx * hx / hy / hy * (lambda_y[j] - 2.0);
    }
    delete[] lambda_y;

    bool is_no_Dirichlet = !(isDirLike(boundary_type_left) || isDirLike(boundary_type_right) ||
                             isDirLike(boundary_type_down) || isDirLike(boundary_type_up));
    bool has_last_vector = true;

    chasing_method_x = new ChasingMethod2D();
    chasing_method_x->init(ny, nx, x_diag, is_no_Dirichlet, has_last_vector, boundary_type_left, boundary_type_right);
}

PoissonSolver2D::~PoissonSolver2D()
{
    delete[] x_diag;
    // delete[] lambda_y;
    delete poisson_fft_y;
    delete chasing_method_x;
}

void PoissonSolver2D::solve(field2& f, bool is_debugmode)
{
    if (env_config && env_config->showCurrentStep)
        std::cout << "[Poisson] solve: start" << std::endl;

    if (env_config && env_config->debugMode && is_debugmode)
    {
        std::string fname_rhs =
            env_config->debugOutputDir + "/rhs_" + domain->name + "_" + std::to_string(solve_call_count);
        IO::field_to_csv(f, fname_rhs);
    }

    buffer.set_size(nx, ny);
    poisson_fft_y->transform(f, buffer);

    buffer.transpose(f);
    buffer.set_size(ny, nx);
    chasing_method_x->chasing(f, buffer);

    buffer.transpose(f);
    buffer.set_size(nx, ny);
    poisson_fft_y->transform_transpose(f, buffer);

    std::swap(f, buffer);

    if (env_config && env_config->debugMode && is_debugmode)
    {
        std::string fname_sol =
            env_config->debugOutputDir + "/sol_" + domain->name + "_" + std::to_string(solve_call_count);
        IO::field_to_csv(f, fname_sol);
    }

    solve_call_count++;
    if (env_config && env_config->showCurrentStep)
    {
        double s_f = f.sum();
        std::cout << "[Poisson] solve: done, f.sum=" << s_f << std::endl;
    }
}

void PoissonSolver2D::cal_lambda() // The current version is only for OpenMP
{
    // This function is calcualting the lambda of the equation with all the diagonal elements equal to 0
    for (int i = 1; i <= ny; i++)
    {
        if (boundary_type_down == PDEBoundaryType::Periodic && boundary_type_up == PDEBoundaryType::Periodic) // P-P
        {
            lambda_y[i - 1] = -2.0 * std::cos(2.0 * pi / ny * std::floor(i / 2.0));
        }
        else if (boundary_type_down == PDEBoundaryType::Neumann && boundary_type_up == PDEBoundaryType::Neumann) // N-N
        {
            lambda_y[i - 1] = -2.0 * std::cos(pi / ny * i);
        }
        else if ((boundary_type_down == PDEBoundaryType::Dirichlet ||
                  boundary_type_down == PDEBoundaryType::Adjacented) &&
                 (boundary_type_up == PDEBoundaryType::Dirichlet ||
                  boundary_type_up == PDEBoundaryType::Adjacented)) // D/Adj - D/Adj
        {
            lambda_y[i - 1] = +2.0 * std::cos(pi / (ny + 1) * i);
        }
        else if (((boundary_type_down == PDEBoundaryType::Dirichlet ||
                   boundary_type_down == PDEBoundaryType::Adjacented) &&
                  boundary_type_up == PDEBoundaryType::Neumann) ||
                 (boundary_type_down == PDEBoundaryType::Neumann &&
                  (boundary_type_up == PDEBoundaryType::Dirichlet ||
                   boundary_type_up == PDEBoundaryType::Adjacented))) // (D/Adj)-N or N-(D/Adj)
        {
            lambda_y[i - 1] = -2.0 * std::cos(2.0 * pi * i / (2 * ny + 1));
        }
    }
}
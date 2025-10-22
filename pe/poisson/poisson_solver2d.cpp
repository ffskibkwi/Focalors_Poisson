#include "poisson_solver2d.h"

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
    : nx(in_domain->nx)
    , ny(in_domain->ny)
    , hx(in_domain->hx)
    , hy(in_domain->hy)
{
    env_config = in_env_config;
    // 新逻辑：必须从 Variable 的 boundary_type_map 读取边界；禁止回退到 Domain
    // if (in_variable == nullptr)
    //     throw std::runtime_error("PoissonSolver2D requires Variable with boundary_type_map; do not use Domain
    //     boundary");

    // auto domIt = in_variable->boundary_type_map.find(in_domain);
    // if (domIt == in_variable->boundary_type_map.end())
    //     throw std::runtime_error("Variable has no boundary map for domain " + in_domain->name);

    // const auto &mp = domIt->second;
    // auto get_required = [&](LocationType loc){
    //     auto it = mp.find(loc);
    //     if (it == mp.end() || it->second == PDEBoundaryType::Null)
    //         throw std::runtime_error("Boundary type missing for domain " + in_domain->name);
    //     return it->second;
    // };

    boundary_type_left  = in_variable->boundary_type_map[in_domain][LocationType::Left];
    boundary_type_right = in_variable->boundary_type_map[in_domain][LocationType::Right];
    boundary_type_down  = in_variable->boundary_type_map[in_domain][LocationType::Down];
    boundary_type_up    = in_variable->boundary_type_map[in_domain][LocationType::Up];

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
    else if ((isDirLike(boundary_type_down) && boundary_type_up == PDEBoundaryType::Neumann) ||
             (boundary_type_down == PDEBoundaryType::Neumann && isDirLike(boundary_type_up)))
    {
        poisson_fft_y = new PoissonFFT2D_DN();
    }
    poisson_fft_y->init(nx, ny);

    lambda_y = new double[ny];
    x_diag   = new double[ny];

    cal_lambda();

    for (int j = 0; j < ny; j++)
    {
        x_diag[j] = -2.0 + hx * hx / hy / hy * lambda_y[j];
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

void PoissonSolver2D::solve(field2& f)
{
    if (env_config && env_config->showCurrentStep)
        std::cout << "[Poisson] solve: start" << std::endl;

    boundary_assembly(f);

    buffer.set_size(nx, ny);
    poisson_fft_y->transform(f, buffer);

    buffer.transpose(f);
    buffer.set_size(ny, nx);
    chasing_method_x->chasing(f, buffer);

    buffer.transpose(f);
    buffer.set_size(nx, ny);
    poisson_fft_y->transform_transpose(f, buffer);

    std::swap(f, buffer);
    if (env_config && env_config->showCurrentStep)
        std::cout << "[Poisson] solve: done" << std::endl;
}

void PoissonSolver2D::cal_lambda() // The current version is only for OpenMP
{
    // This function is calcualting the lambda of the equation with all the diagonal elements equal to 0
    for (int i = 1; i <= ny; i++)
    {
        if (boundary_type_down == PDEBoundaryType::Periodic && boundary_type_up == PDEBoundaryType::Periodic) // P-P
        {
            lambda_y[i - 1] = -2.0 + 2.0 * std::cos(2.0 * pi / ny * std::floor(i / 2.0));
        }
        else if (boundary_type_down == PDEBoundaryType::Neumann && boundary_type_up == PDEBoundaryType::Neumann) // N-N
        {
            lambda_y[i - 1] = -2.0 + 2.0 * std::cos(pi / ny * i);
        }
        else if ((boundary_type_down == PDEBoundaryType::Dirichlet ||
                  boundary_type_down == PDEBoundaryType::Adjacented) &&
                 (boundary_type_up == PDEBoundaryType::Dirichlet ||
                  boundary_type_up == PDEBoundaryType::Adjacented)) // D/Adj - D/Adj
        {
            lambda_y[i - 1] = -2.0 + 2.0 * std::cos(pi / (ny + 1) * i);
        }
        else if (((boundary_type_down == PDEBoundaryType::Dirichlet ||
                   boundary_type_down == PDEBoundaryType::Adjacented) &&
                  boundary_type_up == PDEBoundaryType::Neumann) ||
                 (boundary_type_down == PDEBoundaryType::Neumann &&
                  (boundary_type_up == PDEBoundaryType::Dirichlet ||
                   boundary_type_up == PDEBoundaryType::Adjacented))) // (D/Adj)-N or N-(D/Adj)
        {
            lambda_y[i - 1] = -2.0 + 2.0 * std::cos(2.0 * pi * i / (2 * ny + 1));
        }
    }
}

void PoissonSolver2D::boundary_assembly(field2& f)
{
    auto& var_has_map   = in_variable->has_boundary_value_map[in_domain];
    auto& var_value_map = in_variable->boundary_value_map[in_domain];

    if (boundary_type_left == PDEBoundaryType::Dirichlet && var_has_map[LocationType::Left])
    {
        double* boundary_value = var_value_map[LocationType::Left];
        for (int j = 0; j < ny; j++)
            f(0, j) -= boundary_value[j] / hx / hx;
    }
    if (boundary_type_right == PDEBoundaryType::Dirichlet && var_has_map[LocationType::Right])
    {
        double* boundary_value = var_value_map[LocationType::Right];
        for (int j = 0; j < ny; j++)
            f(nx - 1, j) -= boundary_value[j] / hx / hx;
    }

    if (boundary_type_down == PDEBoundaryType::Dirichlet && var_has_map[LocationType::Down])
    {
        double* boundary_value = var_value_map[LocationType::Down];
        for (int i = 0; i < nx; i++)
            f(i, 0) -= boundary_value[i] / hy / hy;
    }
    if (boundary_type_up == PDEBoundaryType::Dirichlet && var_has_map[LocationType::Up])
    {
        double* boundary_value = var_value_map[LocationType::Up];
        for (int i = 0; i < nx; i++)
            f(i, ny - 1) -= boundary_value[i] / hy / hy;
    }

    if (boundary_type_left == PDEBoundaryType::Neumann && var_has_map[LocationType::Left])
    {
        double* boundary_value = var_value_map[LocationType::Left];
        for (int j = 0; j < ny; j++)
            f(0, j) += boundary_value[j] / hx;
    }
    if (boundary_type_right == PDEBoundaryType::Neumann && var_has_map[LocationType::Right])
    {
        double* boundary_value = var_value_map[LocationType::Left];
        for (int j = 0; j < ny; j++)
            f(nx - 1, j) -= boundary_value[j] / hx;
    }
    if (boundary_type_down == PDEBoundaryType::Neumann && var_has_map[LocationType::Down])
    {
        double* boundary_value = var_value_map[LocationType::Down];
        for (int i = 0; i < nx; i++)
            f(i, 0) += boundary_value[i] / hy;
    }
    if (boundary_type_up == PDEBoundaryType::Neumann && var_has_map[LocationType::Up])
    {
        double* boundary_value = var_value_map[LocationType::Up];
        for (int i = 0; i < nx; i++)
            f(i, ny - 1) -= boundary_value[i] / hy;
    }
}
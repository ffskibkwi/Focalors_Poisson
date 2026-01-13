#include "poisson_solver3d.h"

PoissonSolver3D::PoissonSolver3D(int                in_nx,
                                 int                in_ny,
                                 int                in_nz,
                                 double             in_hx,
                                 double             in_hy,
                                 double             in_hz,
                                 PDEBoundaryType    in_boundary_type_left,
                                 PDEBoundaryType    in_boundary_type_right,
                                 PDEBoundaryType    in_boundary_type_front,
                                 PDEBoundaryType    in_boundary_type_back,
                                 PDEBoundaryType    in_boundary_type_down,
                                 PDEBoundaryType    in_boundary_type_up,
                                 EnvironmentConfig* in_env_config)
    : nx(in_nx)
    , ny(in_ny)
    , nz(in_nz)
    , hx(in_hx)
    , hy(in_hy)
    , hz(in_hz)
    , boundary_type_left(in_boundary_type_left)
    , boundary_type_right(in_boundary_type_right)
    , boundary_type_front(in_boundary_type_front)
    , boundary_type_back(in_boundary_type_back)
    , boundary_type_down(in_boundary_type_down)
    , boundary_type_up(in_boundary_type_up)
{
    env_config = in_env_config;
    init();
}

PoissonSolver3D::PoissonSolver3D(Domain3DUniform* in_domain, Variable* in_variable, EnvironmentConfig* in_env_config)
    : domain(in_domain)
    , var(in_variable)
    , nx(in_domain->nx)
    , ny(in_domain->ny)
    , nz(in_domain->nz)
    , hx(in_domain->hx)
    , hy(in_domain->hy)
    , hz(in_domain->hz)
{
    env_config = in_env_config;
    // 新逻辑：必须从 Variable 的 boundary_type_map 读取边界；禁止回退到 Domain
    // if (in_variable == nullptr)
    //     throw std::runtime_error("PoissonSolver3D requires Variable with boundary_type_map; do not use Domain
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

    boundary_type_left  = var->boundary_type_map[domain][LocationType::Left];
    boundary_type_right = var->boundary_type_map[domain][LocationType::Right];
    boundary_type_front = var->boundary_type_map[domain][LocationType::Front];
    boundary_type_back  = var->boundary_type_map[domain][LocationType::Back];
    boundary_type_down  = var->boundary_type_map[domain][LocationType::Down];
    boundary_type_up    = var->boundary_type_map[domain][LocationType::Up];

    init();
}

void PoissonSolver3D::init()
{
    buffer.init(nx,
                ny,
                nz,
                "buffer"); // Here init a new field for buffer, but it can be replaced by another existed buffer field

    auto isDirLike = [](PDEBoundaryType t) {
        return t == PDEBoundaryType::Dirichlet || t == PDEBoundaryType::Adjacented;
    };

    auto create_fft = [&](PoissonFFT3D* fft, PDEBoundaryType bound_neg, PDEBoundaryType bound_pos) {
        if (bound_neg == PDEBoundaryType::Periodic && bound_pos == PDEBoundaryType::Periodic)
        {
            fft = new PoissonFFT3D_PP();
        }
        else if (bound_neg == PDEBoundaryType::Neumann && bound_pos == PDEBoundaryType::Neumann)
        {
            fft = new PoissonFFT3D_NN();
        }
        else if (isDirLike(bound_neg) && isDirLike(bound_pos))
        {
            fft = new PoissonFFT3D_DD();
        }
        else if (isDirLike(bound_neg) && bound_pos == PDEBoundaryType::Neumann)
        {
            fft = new PoissonFFT3D_DN();
        }
        else if (bound_neg == PDEBoundaryType::Neumann && isDirLike(bound_pos))
        {
            fft = new PoissonFFT3D_ND();
        }
        fft->init(nx, ny, nz);
    };
    create_fft(poisson_fft_z, boundary_type_down, boundary_type_up);
    create_fft(poisson_fft_y, boundary_type_front, boundary_type_back);

    double* lambda_z = new double[nz];
    double* lambda_y = new double[ny];
    x_diag           = new double*[nz];

    cal_lambda(lambda_z, nz, 1, nz, boundary_type_down, boundary_type_up);
    cal_lambda(lambda_y, ny, 1, ny, boundary_type_front, boundary_type_back);

    for (int k = 0; k < nz; k++)
    {
        x_diag[k] = new double[ny];
        for (int j = 0; j < ny; j++)
        {
            x_diag[k][j] = -2.0 + hx * hx / hy / hy * lambda_y[j] + hx * hx / hz / hz * lambda_z[k];
        }
    }
    delete[] lambda_y;
    delete[] lambda_z;

    bool is_no_Dirichlet =
        !(isDirLike(boundary_type_left) || isDirLike(boundary_type_right) || isDirLike(boundary_type_front) ||
          isDirLike(boundary_type_back) || isDirLike(boundary_type_down) || isDirLike(boundary_type_up));
    bool has_last_vector = true;

    chasing_method_x = new ChasingMethod3D();
    chasing_method_x->init(
        nz, ny, nx, x_diag, is_no_Dirichlet, has_last_vector, boundary_type_left, boundary_type_right);
}

PoissonSolver3D::~PoissonSolver3D()
{
    for (int k = 0; k < nz; k++)
    {
        delete[] x_diag[k];
    }
    delete[] x_diag;
    delete poisson_fft_z;
    delete poisson_fft_y;
    delete chasing_method_x;
}

void PoissonSolver3D::solve(field3& f)
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
    {
        double s_f = f.sum();
        std::cout << "[Poisson] solve: done, f.sum=" << s_f << std::endl;
    }
}

void PoissonSolver3D::cal_lambda(double*           lambda,
                                 int               global_length,
                                 int               begin,
                                 int               local_length,
                                 FluidBoundaryType BoundTypeNegative,
                                 FluidBoundaryType BoundTypePositive) // The current version is only for OpenMP
{
    // This function is calcualting the lambda of the equation with all the diagonal elements equal to 0
    for (int i = begin; i < begin + local_length; i++)
    {
        if (boundary_type_down == PDEBoundaryType::Periodic && boundary_type_up == PDEBoundaryType::Periodic) // P-P
        {
            lambda[i - begin] = -2.0 - 2.0 * std::cos(2.0 * pi / global_length * std::floor(i / 2.0));
        }
        else if (boundary_type_down == PDEBoundaryType::Neumann && boundary_type_up == PDEBoundaryType::Neumann) // N-N
        {
            lambda[i - begin] = -2.0 - 2.0 * std::cos(pi / global_length * i);
        }
        else if ((boundary_type_down == PDEBoundaryType::Dirichlet ||
                  boundary_type_down == PDEBoundaryType::Adjacented) &&
                 (boundary_type_up == PDEBoundaryType::Dirichlet ||
                  boundary_type_up == PDEBoundaryType::Adjacented)) // D/Adj - D/Adj
        {
            lambda[i - begin] = -2.0 + 2.0 * std::cos(pi / (global_length + 1) * i);
        }
        else if (((boundary_type_down == PDEBoundaryType::Dirichlet ||
                   boundary_type_down == PDEBoundaryType::Adjacented) &&
                  boundary_type_up == PDEBoundaryType::Neumann) ||
                 (boundary_type_down == PDEBoundaryType::Neumann &&
                  (boundary_type_up == PDEBoundaryType::Dirichlet ||
                   boundary_type_up == PDEBoundaryType::Adjacented))) // (D/Adj)-N or N-(D/Adj)
        {
            lambda[i - begin] = -2.0 - 2.0 * std::cos(2.0 * pi * i / (2 * global_length + 1));
        }
    }
}

void PoissonSolver3D::boundary_assembly(field3& f)
{
    auto& var_has_map   = var->has_boundary_value_map[domain];
    auto& var_value_map = var->boundary_value_map[domain];

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
        double* boundary_value = var_value_map[LocationType::Right];
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
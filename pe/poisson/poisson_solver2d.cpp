#include "poisson_solver2d.h"

PoissonSolver2D::PoissonSolver2D(int in_nx, int in_ny, double in_hx, double in_hy, PDEBoundaryType in_BoundaryTypeXNegative, PDEBoundaryType in_BoundaryTypeXPositive, PDEBoundaryType in_BoundaryTypeYNegative, PDEBoundaryType in_BoundaryTypeYPositive, EnvironmentConfig* in_env_config)
    : nx(in_nx), 
    ny(in_ny), 
    hx(in_hx), 
    hy(in_hy), 
    BoundaryTypeXNegative(in_BoundaryTypeXNegative), 
    BoundaryTypeXPositive(in_BoundaryTypeXPositive), 
    BoundaryTypeYNegative(in_BoundaryTypeYNegative), 
    BoundaryTypeYPositive(in_BoundaryTypeYPositive)
{
    env_config = in_env_config;
    init();
}

PoissonSolver2D::PoissonSolver2D(Domain2DUniform* in_domain, Variable* in_variable, EnvironmentConfig* in_env_config)
    : nx(in_domain->nx), 
    ny(in_domain->ny), 
    hx(in_domain->hx), 
    hy(in_domain->hy)
{
    env_config = in_env_config;
    // 新逻辑：必须从 Variable 的 boundary_type_map 读取边界；禁止回退到 Domain
    if (in_variable == nullptr)
        throw std::runtime_error("PoissonSolver2D requires Variable with boundary_type_map; do not use Domain boundary");

    auto domIt = in_variable->boundary_type_map.find(in_domain);
    if (domIt == in_variable->boundary_type_map.end())
        throw std::runtime_error("Variable has no boundary map for domain " + in_domain->name);

    const auto &mp = domIt->second;
    auto get_required = [&](LocationType loc){
        auto it = mp.find(loc);
        if (it == mp.end() || it->second == PDEBoundaryType::Null)
            throw std::runtime_error("Boundary type missing for domain " + in_domain->name);
        return it->second;
    };

    BoundaryTypeXNegative = get_required(LocationType::Left);
    BoundaryTypeXPositive = get_required(LocationType::Right);
    BoundaryTypeYNegative = get_required(LocationType::Down);
    BoundaryTypeYPositive = get_required(LocationType::Up);

    init();
}

void PoissonSolver2D::init()
{
    buffer.init(nx, ny, "buffer");    //Here init a new field for buffer, but it can be replaced by another existed buffer field

    auto isDirLike = [](PDEBoundaryType t){ return t == PDEBoundaryType::Dirichlet || t == PDEBoundaryType::Adjacented; };

    if (BoundaryTypeYNegative == PDEBoundaryType::Periodic && BoundaryTypeYPositive == PDEBoundaryType::Periodic)
    {
        poisson_fft_y = new PoissonFFT2D_PP();
    }
    else if (BoundaryTypeYNegative == PDEBoundaryType::Neumann && BoundaryTypeYPositive == PDEBoundaryType::Neumann)
    {
        poisson_fft_y = new PoissonFFT2D_NN();
    }
    else if (isDirLike(BoundaryTypeYNegative) && isDirLike(BoundaryTypeYPositive))
    {
        poisson_fft_y = new PoissonFFT2D_DD();
    }
    else if ((isDirLike(BoundaryTypeYNegative) && BoundaryTypeYPositive == PDEBoundaryType::Neumann) ||
        (BoundaryTypeYNegative == PDEBoundaryType::Neumann && isDirLike(BoundaryTypeYPositive)))
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

    bool is_no_Dirichlet = !(isDirLike(BoundaryTypeXNegative) || isDirLike(BoundaryTypeXPositive) ||
            isDirLike(BoundaryTypeYNegative) || isDirLike(BoundaryTypeYPositive));
    bool has_last_vector = true;

    chasing_method_x = new ChasingMethod2D();
    chasing_method_x->init(ny, nx, x_diag, is_no_Dirichlet, has_last_vector, BoundaryTypeXNegative, BoundaryTypeXPositive);
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

void PoissonSolver2D::cal_lambda()  //The current version is only for OpenMP
{
    // This function is calcualting the lambda of the equation with all the diagonal elements equal to 0
    for (int i = 1; i <= ny; i++)
    {
        if (BoundaryTypeYNegative == PDEBoundaryType::Periodic &&
        BoundaryTypeYPositive == PDEBoundaryType::Periodic) // P-P
        {
        lambda_y[i - 1] = -2.0 + 2.0 * std::cos(2.0 * pi / ny * std::floor(i / 2.0));
        }
        else if (BoundaryTypeYNegative == PDEBoundaryType::Neumann &&
            BoundaryTypeYPositive == PDEBoundaryType::Neumann) // N-N
        {
        lambda_y[i - 1] = -2.0 + 2.0 * std::cos(pi / ny * i);
        }
        else if ((BoundaryTypeYNegative == PDEBoundaryType::Dirichlet || BoundaryTypeYNegative == PDEBoundaryType::Adjacented) &&
            (BoundaryTypeYPositive == PDEBoundaryType::Dirichlet || BoundaryTypeYPositive == PDEBoundaryType::Adjacented)) // D/Adj - D/Adj
        {
        lambda_y[i - 1] = -2.0 + 2.0 * std::cos(pi / (ny + 1) * i);
        }
        else if (((BoundaryTypeYNegative == PDEBoundaryType::Dirichlet || BoundaryTypeYNegative == PDEBoundaryType::Adjacented) && BoundaryTypeYPositive == PDEBoundaryType::Neumann) ||
            (BoundaryTypeYNegative == PDEBoundaryType::Neumann && (BoundaryTypeYPositive == PDEBoundaryType::Dirichlet || BoundaryTypeYPositive == PDEBoundaryType::Adjacented))) // (D/Adj)-N or N-(D/Adj)
        {
        lambda_y[i - 1] = -2.0 + 2.0 * std::cos(2.0 * pi * i / (2 * ny + 1));
        }
    }
}
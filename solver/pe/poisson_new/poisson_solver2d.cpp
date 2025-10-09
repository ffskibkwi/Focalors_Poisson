#include "poisson_solver2d.h"

PoissonSolver2D::PoissonSolver2D(int in_nx, int in_ny, double in_hx, double in_hy, PDEBoundaryType in_BoundaryTypeXNegative, PDEBoundaryType in_BoundaryTypeXPositive, PDEBoundaryType in_BoundaryTypeYNegative, PDEBoundaryType in_BoundaryTypeYPositive)
    : nx(in_nx), 
    ny(in_ny), 
    hx(in_hx), 
    hy(in_hy), 
    BoundaryTypeXNegative(in_BoundaryTypeXNegative), 
    BoundaryTypeXPositive(in_BoundaryTypeXPositive), 
    BoundaryTypeYNegative(in_BoundaryTypeYNegative), 
    BoundaryTypeYPositive(in_BoundaryTypeYPositive)
{
    init();
}

PoissonSolver2D::PoissonSolver2D(Domain2DUniform& in_domain)
    : nx(in_domain.nx), 
    ny(in_domain.ny), 
    hx(in_domain.hx), 
    hy(in_domain.hy), 
    BoundaryTypeXNegative(in_domain.BoundaryTypeXNegative), 
    BoundaryTypeXPositive(in_domain.BoundaryTypeXPositive), 
    BoundaryTypeYNegative(in_domain.BoundaryTypeYNegative), 
    BoundaryTypeYPositive(in_domain.BoundaryTypeYPositive)
{
    init();
}

void PoissonSolver2D::init()
{
    buffer.init(nx, ny, "buffer");    //Here init a new field for buffer, but it can be replaced by another existed buffer field

    if (BoundaryTypeYNegative == PDEBoundaryType::Periodic && BoundaryTypeYPositive == PDEBoundaryType::Periodic)
    {
        poisson_fft_y = new PoissonFFT2D_PP();
    }
    else if (BoundaryTypeYNegative == PDEBoundaryType::Neumann && BoundaryTypeYPositive == PDEBoundaryType::Neumann)
    {
        poisson_fft_y = new PoissonFFT2D_NN();
    }
    else if (BoundaryTypeYNegative == PDEBoundaryType::Dirichlet && BoundaryTypeYPositive == PDEBoundaryType::Dirichlet)
    {
        poisson_fft_y = new PoissonFFT2D_DD();
    }
    else if ((BoundaryTypeYNegative == PDEBoundaryType::Dirichlet && BoundaryTypeYPositive == PDEBoundaryType::Neumann) ||
        (BoundaryTypeYNegative == PDEBoundaryType::Neumann && BoundaryTypeYPositive == PDEBoundaryType::Dirichlet))
    {
        poisson_fft_y = new PoissonFFT2D_DN();
    }
    poisson_fft_y->init(nx, ny);


    lambda_y = new double[ny];
    x_diag   = new double[ny];

    cal_lambda();
    // cal_lambda(lambda_y, ny, 1, ny, BoundaryTypeYNegative, BoundaryTypeYPositive);

    for (int j = 0; j < ny; j++)
    {
        x_diag[j] = -2.0 + hx * hx / hy / hy * (lambda_y[j] - 2.0);
    }
    delete[] lambda_y;

    bool is_no_Dirichlet =
        (BoundaryTypeXNegative != PDEBoundaryType::Dirichlet && BoundaryTypeXPositive != PDEBoundaryType::Dirichlet &&
            BoundaryTypeYNegative != PDEBoundaryType::Dirichlet && BoundaryTypeYPositive != PDEBoundaryType::Dirichlet);
    bool has_last_vector = true;

    chasing_method_x = new ChasingMethod2D();
    chasing_method_x->init(nx, ny, x_diag, is_no_Dirichlet, has_last_vector, BoundaryTypeXNegative, BoundaryTypeXPositive);
}

void PoissonSolver2D::~PoissonSolver2D()
{
    delete[] x_diag;
    // delete[] lambda_y;
    delete poisson_fft_y;
    delete chasing_method_x;
}

void PoissonSolver2D::solve(field2& f)
{
    buffer.set_size(nx, ny);
    poisson_fft_y->transform(f, buffer);

    chasing_method_x->chasing(buffer, f);

    buffer.set_size(nx, ny);
    poisson_fft_y->transform_transpose(f, buffer);

    std::swap(f, buffer);
}

void PoissonSolver2D::cal_lambda()  //The current version is only for OpenMP
{
    // This function is calcualting the lambda of the equation with all the diagonal elements equal to 0
    for (int i = 0; i < ny; i++)
    {
        if (BoundaryTypeYNegative == PDEBoundaryType::Periodic &&
        BoundaryTypeYPositive == PDEBoundaryType::Periodic) // P-P
        {
        lambda[i] = -2.0 * std::cos(2.0 * pi / ny * std::floor(i / 2.0));
        }
        else if (BoundaryTypeYNegative == PDEBoundaryType::Neumann &&
            BoundaryTypeYPositive == PDEBoundaryType::Neumann) // N-N
        {
        lambda[i] = -2.0 * std::cos(pi / ny * i);
        }
        else if (BoundaryTypeYNegative == PDEBoundaryType::Dirichlet &&
            BoundaryTypeYPositive == PDEBoundaryType::Dirichlet) // D-D
        {
        lambda[i] = -2.0 * std::cos(pi / (ny + 1) * i);
        }
        else if ((BoundaryTypeYNegative == PDEBoundaryType::Dirichlet && BoundaryTypeYPositive == PDEBoundaryType::Neumann) ||
            (BoundaryTypeYNegative == PDEBoundaryType::Neumann && BoundaryTypeYPositive == PDEBoundaryType::Dirichlet)) // D-N or N-D
        {
        lambda[i] = -2.0 * std::cos(2.0 * pi * i / (2 * ny + 1));
        }
    }
}
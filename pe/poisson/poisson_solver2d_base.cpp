#include "poisson_solver2d_base.h"

PoissonSolver2DBase::PoissonSolver2DBase(int             in_nx,
                                         int             in_ny,
                                         double          in_hx,
                                         double          in_hy,
                                         PDEBoundaryType in_boundary_type_xneg,
                                         PDEBoundaryType in_boundary_type_xpos,
                                         PDEBoundaryType in_boundary_type_yneg,
                                         PDEBoundaryType in_boundary_type_ypos)
    : nx(in_nx)
    , ny(in_ny)
    , hx(in_hx)
    , hy(in_hy)
    , boundary_type_xneg(in_boundary_type_xneg)
    , boundary_type_xpos(in_boundary_type_xpos)
    , boundary_type_yneg(in_boundary_type_yneg)
    , boundary_type_ypos(in_boundary_type_ypos)
{}

PoissonSolver2DBase::~PoissonSolver2DBase()
{
    delete[] x_diag;
    delete poisson_fft_y;
    delete chasing_method_x;
}

void PoissonSolver2DBase::cal_lambda(double*         lambda,
                                     int             global_length,
                                     int             begin,
                                     int             local_length,
                                     PDEBoundaryType BoundTypeNegative,
                                     PDEBoundaryType BoundTypePositive)
{
    for (int i = begin; i < begin + local_length; i++)
    {
        if (BoundTypeNegative == PDEBoundaryType::Periodic && BoundTypePositive == PDEBoundaryType::Periodic) // P-P
        {
            lambda[i - begin] = -2.0 * std::cos(2.0 * pi / global_length * std::floor(i / 2.0));
        }
        else if (BoundTypeNegative == PDEBoundaryType::Neumann && BoundTypePositive == PDEBoundaryType::Neumann) // N-N
        {
            lambda[i - begin] = -2.0 * std::cos(pi / global_length * i);
        }
        else if ((BoundTypeNegative == PDEBoundaryType::Dirichlet ||
                  BoundTypeNegative == PDEBoundaryType::Adjacented) &&
                 (BoundTypePositive == PDEBoundaryType::Dirichlet ||
                  BoundTypePositive == PDEBoundaryType::Adjacented)) // D/Adj - D/Adj
        {
            lambda[i - begin] = +2.0 * std::cos(pi / (global_length + 1) * i);
        }
        else if (((BoundTypeNegative == PDEBoundaryType::Dirichlet ||
                   BoundTypeNegative == PDEBoundaryType::Adjacented) &&
                  BoundTypePositive == PDEBoundaryType::Neumann) ||
                 (BoundTypeNegative == PDEBoundaryType::Neumann &&
                  (BoundTypePositive == PDEBoundaryType::Dirichlet ||
                   BoundTypePositive == PDEBoundaryType::Adjacented))) // (D/Adj)-N or N-(D/Adj)
        {
            lambda[i - begin] = -2.0 * std::cos(2.0 * pi * i / (2 * global_length + 1));
        }
    }
}

void PoissonSolver2DBase::create_fft(PoissonFFT2D*&  fft,
                                     PDEBoundaryType bound_neg,
                                     PDEBoundaryType bound_pos,
                                     int             n1,
                                     int             n2)
{
    if (bound_neg == PDEBoundaryType::Periodic && bound_pos == PDEBoundaryType::Periodic)
    {
        fft = new PoissonFFT2D_PP();
    }
    else if (bound_neg == PDEBoundaryType::Neumann && bound_pos == PDEBoundaryType::Neumann)
    {
        fft = new PoissonFFT2D_NN();
    }
    else if (isDirLike(bound_neg) && isDirLike(bound_pos))
    {
        fft = new PoissonFFT2D_DD();
    }
    else if (isDirLike(bound_neg) && bound_pos == PDEBoundaryType::Neumann)
    {
        fft = new PoissonFFT2D_DN();
    }
    else if (bound_neg == PDEBoundaryType::Neumann && isDirLike(bound_pos))
    {
        fft = new PoissonFFT2D_ND();
    }
    fft->init(n1, n2);
}
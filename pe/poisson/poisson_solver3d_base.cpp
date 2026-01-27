#include "poisson_solver3d_base.h"

PoissonSolver3DBase::PoissonSolver3DBase(int             in_nx,
                                         int             in_ny,
                                         int             in_nz,
                                         double          in_hx,
                                         double          in_hy,
                                         double          in_hz,
                                         PDEBoundaryType in_boundary_type_left,
                                         PDEBoundaryType in_boundary_type_right,
                                         PDEBoundaryType in_boundary_type_front,
                                         PDEBoundaryType in_boundary_type_back,
                                         PDEBoundaryType in_boundary_type_down,
                                         PDEBoundaryType in_boundary_type_up)
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
{}

PoissonSolver3DBase::~PoissonSolver3DBase()
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

void PoissonSolver3DBase::cal_lambda(double*         lambda,
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
            lambda[i - begin] = -2.0 - 2.0 * std::cos(2.0 * pi / global_length * std::floor(i / 2.0));
        }
        else if (BoundTypeNegative == PDEBoundaryType::Neumann && BoundTypePositive == PDEBoundaryType::Neumann) // N-N
        {
            lambda[i - begin] = -2.0 - 2.0 * std::cos(pi / global_length * i);
        }
        else if ((BoundTypeNegative == PDEBoundaryType::Dirichlet ||
                  BoundTypeNegative == PDEBoundaryType::Adjacented) &&
                 (BoundTypePositive == PDEBoundaryType::Dirichlet ||
                  BoundTypePositive == PDEBoundaryType::Adjacented)) // D/Adj - D/Adj
        {
            lambda[i - begin] = -2.0 + 2.0 * std::cos(pi / (global_length + 1) * i);
        }
        else if (((BoundTypeNegative == PDEBoundaryType::Dirichlet ||
                   BoundTypeNegative == PDEBoundaryType::Adjacented) &&
                  BoundTypePositive == PDEBoundaryType::Neumann) ||
                 (BoundTypeNegative == PDEBoundaryType::Neumann &&
                  (BoundTypePositive == PDEBoundaryType::Dirichlet ||
                   BoundTypePositive == PDEBoundaryType::Adjacented))) // (D/Adj)-N or N-(D/Adj)
        {
            lambda[i - begin] = -2.0 - 2.0 * std::cos(2.0 * pi * i / (2 * global_length + 1));
        }
    }
}

void PoissonSolver3DBase::create_fft(PoissonFFT3D*&  fft,
                                     PDEBoundaryType bound_neg,
                                     PDEBoundaryType bound_pos,
                                     int             n1,
                                     int             n2,
                                     int             n3)
{
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
    fft->init(n1, n2, n3);
}
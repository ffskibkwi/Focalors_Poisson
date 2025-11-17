#include "chasing_method2d.h"

void ChasingMethod2D::init(int nx_in, int ny_in, double* _x_diag, bool _is_no_Dirichlet, bool _has_last_vector,
                           PDEBoundaryType boundary_type_down, PDEBoundaryType boundary_type_up)
{
    nx               = nx_in;
    ny               = ny_in;
    x_diag           = _x_diag;
    is_no_Dirichlet  = _is_no_Dirichlet;
    has_last_vector  = _has_last_vector;
    this->boundary_type_down = boundary_type_down;
    this->boundary_type_up = boundary_type_up;

    y.init(nx, ny, "y");

    if (boundary_type_down == PDEBoundaryType::Periodic && boundary_type_up == PDEBoundaryType::Periodic)
    {
        z_1.init(nx, ny, "z_1");
        z_2.init(nx, ny, "z_2");
        c_diri.init(nx, ny - 1, "c_diri");
        y_diri.init(nx, ny, "y_diri");
        for (int i = 0; i < nx; i++)
        {
            chasing_periodic_precompute(x_diag[i], ny, z_1.get_ptr(i), z_2.get_ptr(i), c_diri.get_ptr(i));
        }
    }
    else
    {
        c.init(nx, ny - 1, "c");
        for (int i = 0; i < nx; i++)
        {
            chasing_standard_precompute(x_diag[i], ny, c.get_ptr(i), boundary_type_down, boundary_type_up);
        }
    }
}

void ChasingMethod2D::chasing(field2& f, field2& p)
{
    if (boundary_type_down == PDEBoundaryType::Periodic && boundary_type_up == PDEBoundaryType::Periodic)
    {
        for (int i = 0; i < nx; i++)
        {
            if (is_no_Dirichlet && has_last_vector && i == (nx - 1))
            {
                chasing_periodic_singular(ny, f.get_ptr(i, 0), p.get_ptr(i, 0));
            }
            else
            {
                chasing_periodic_with_precompute(x_diag[i],
                                                 ny,
                                                 f.get_ptr(i),
                                                 c_diri.get_ptr(i),
                                                 y_diri.get_ptr(i),
                                                 z_1.get_ptr(i),
                                                 z_2.get_ptr(i),
                                                 y.get_ptr(i),
                                                 p.get_ptr(i));
            }
        }
    }
    else
    {
        for (int i = 0; i < nx; i++)
        {
            if (is_no_Dirichlet && has_last_vector && i == (nx - 1))
            {
                chasing_standard_singular(ny, f.get_ptr(i), p.get_ptr(i));
            }
            else
            {
                chasing_standard_without_precompute(x_diag[i],
                                                    ny,
                                                    f.get_ptr(i),
                                                    c.get_ptr(i),
                                                    y.get_ptr(i),
                                                    p.get_ptr(i),
                                                    boundary_type_down,
                                                    boundary_type_up);
            }
        }
    }
}



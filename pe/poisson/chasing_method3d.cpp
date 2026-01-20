#include "chasing_method3d.h"

void ChasingMethod3D::init(int             nx_in,
                           int             ny_in,
                           int             nz_in,
                           double**        _x_diag,
                           bool            _is_no_Dirichlet,
                           bool            _has_last_vector,
                           PDEBoundaryType boundary_type_down,
                           PDEBoundaryType boundary_type_up)
{
    nx                       = nx_in;
    ny                       = ny_in;
    nz                       = nz_in;
    x_diag                   = _x_diag;
    is_no_Dirichlet          = _is_no_Dirichlet;
    has_last_vector          = _has_last_vector;
    this->boundary_type_down = boundary_type_down;
    this->boundary_type_up   = boundary_type_up;

    y.init(nx, ny, nz, "y");

    if (boundary_type_down == PDEBoundaryType::Periodic && boundary_type_up == PDEBoundaryType::Periodic)
    {
        z_1.init(nx, ny, nz, "z_1");
        z_2.init(nx, ny, nz, "z_2");
        c_diri.init(nx, ny, nz - 1, "c_diri");
        y_diri.init(nx, ny, nz, "y_diri");
        for (int i = 0; i < nx; i++)
        {
            for (int j = 0; j < ny; j++)
            {
                chasing_periodic_precompute(
                    x_diag[i][j], nz, z_1.get_ptr(i, j, 0), z_2.get_ptr(i, j, 0), c_diri.get_ptr(i, j, 0));
            }
        }
    }
    else
    {
        c.init(nx, ny, nz - 1, "c");
        for (int i = 0; i < nx; i++)
        {
            for (int j = 0; j < ny; j++)
            {
                chasing_standard_precompute(x_diag[i][j], nz, c.get_ptr(i, j, 0), boundary_type_down, boundary_type_up);
            }
        }
    }
}

void ChasingMethod3D::chasing(field3& f, field3& p)
{
    if (boundary_type_down == PDEBoundaryType::Periodic && boundary_type_up == PDEBoundaryType::Periodic)
    {
        for (int i = 0; i < nx; i++)
        {
            for (int j = 0; j < ny; j++)
            {
                if (is_no_Dirichlet && has_last_vector && i == (nx - 1))
                {
                    chasing_periodic_singular(nz, f.get_ptr(i, j, 0), p.get_ptr(i, j, 0));
                }
                else
                {
                    chasing_periodic_with_precompute(x_diag[i][j],
                                                     nz,
                                                     f.get_ptr(i, j, 0),
                                                     c_diri.get_ptr(i, j, 0),
                                                     y_diri.get_ptr(i, j, 0),
                                                     z_1.get_ptr(i, j, 0),
                                                     z_2.get_ptr(i, j, 0),
                                                     y.get_ptr(i, j, 0),
                                                     p.get_ptr(i, j, 0));
                }
            }
        }
    }
    else
    {
        for (int i = 0; i < nx; i++)
        {
            for (int j = 0; j < ny; j++)
            {
                if (is_no_Dirichlet && has_last_vector && i == (nx - 1) && j == (ny - 1))
                {
                    chasing_standard_singular(nz, f.get_ptr(i, j, 0), p.get_ptr(i, j, 0));
                }
                else
                {
                    chasing_standard_without_precompute(x_diag[i][j],
                                                        nz,
                                                        f.get_ptr(i, j, 0),
                                                        c.get_ptr(i, j, 0),
                                                        y.get_ptr(i, j, 0),
                                                        p.get_ptr(i, j, 0),
                                                        boundary_type_down,
                                                        boundary_type_up);
                }
            }
        }
    }
}

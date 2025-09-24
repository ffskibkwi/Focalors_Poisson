#pragma once

#include "pch.h"

#include "core/boundary/boundary_type.h"
#include "core/parallel/mpi/mpi_misc.h"

#include "chasing_method_base.h"

/**
 * @brief chasing method in z direction on xyz field.
 *
 * @tparam BoundTypeZNegative fluid boundary type in z negative.
 * @tparam BoundTypeZPositive fluid boundary type in z positive.
 * @note
 * In this class, the chasing method is defaultly considered that implement on x-direction, i.e.
 * the chasing function implenent on A(i,j)p(i,j,.)=f(i,j,.) for each i in (1, Nx) and j in (1, Ny).
 * Therefore, in the FFT-based algorithm, there is a conversion from xyz to zyx.
 * But I think it is easy to comprehence. ^-^
 */
template<FluidBoundaryType BoundTypeZNegative, FluidBoundaryType BoundTypeZPositive>
class ChasingMethod3D : public ChasingMethodBase
{
public:
    void init(int nx_in, int ny_in, int nz_in, double** _x_diag, bool _is_no_Dirichlet, bool _has_last_vector)
    {
        nx              = nx_in;
        ny              = ny_in;
        nz              = nz_in;
        x_diag          = _x_diag;
        is_no_Dirichlet = _is_no_Dirichlet;
        has_last_vector = _has_last_vector;

        y.init(nx, ny, nz, "y");

        if constexpr (BoundTypeZNegative == FluidBoundaryType::Periodic &&
                      BoundTypeZPositive == FluidBoundaryType::Periodic)
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
                    chasing_standard_precompute(
                        x_diag[i][j], nz, c.get_ptr(i, j, 0), BoundTypeZNegative, BoundTypeZPositive);
                }
            }
        }
    }

    void chasing(field3& f, field3& p)
    {
        if constexpr (BoundTypeZNegative == FluidBoundaryType::Periodic &&
                      BoundTypeZPositive == FluidBoundaryType::Periodic)
        {
            for (int i = 0; i < nx; i++)
            {
                for (int j = 0; j < ny; j++)
                {
                    if (is_no_Dirichlet && has_last_vector && i == (nx - 1) && j == (ny - 1))
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
                                                            BoundTypeZNegative,
                                                            BoundTypeZPositive);
                    }
                }
            }
        }
    }

private:
    double** x_diag          = nullptr;
    bool     is_no_Dirichlet = false;
    bool     has_last_vector = true;
    int      nx, ny, nz;

    // Intermediate vectors
    field3 y;                        // For all boundary types
    field3 z_1, z_2, c_diri, y_diri; // For periodic boundary types
    field3 c;                        // For standard boundary types
};

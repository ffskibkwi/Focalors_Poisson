#pragma once

#include "base/pch.h"

#include "base/location_boundary.h"
#include "base/parallel/mpi/mpi_misc.h"

#include "chasing_method_base.h"

/**
 * @brief chasing method in y direction on xy field.
 *
 * @tparam BoundTypeYNegative fluid boundary type in y negative.
 * @tparam BoundTypeYPositive fluid boundary type in y positive.
 * @note
 * In this class, the chasing method is defaultly considered that implement on x-direction, i.e.
 * the chasing function implenent on A(i,j)p(i,j,.)=f(i,j,.) for each i in (1, Nx) and j in (1, Ny).
 * Therefore, in the FFT-based algorithm, there is a conversion from xyz to zyx.
 * But I think it is easy to comprehence. ^-^
 */
template<PDEBoundaryType BoundTypeYNegative, PDEBoundaryType BoundTypeYPositive>
class ChasingMethod2D : public ChasingMethodBase
{
public:
    void init(int nx_in, int ny_in, double* _x_diag, bool _is_no_Dirichlet, bool _has_last_vector)
    {
        nx              = nx_in;
        ny              = ny_in;
        x_diag          = _x_diag;
        is_no_Dirichlet = _is_no_Dirichlet;
        has_last_vector = _has_last_vector;

        y.init(nx, ny, "y");

        if constexpr (BoundTypeYNegative == PDEBoundaryType::Periodic &&
                      BoundTypeYPositive == PDEBoundaryType::Periodic)
        {
            z_1.init(nx, ny, "z_1");
            z_2.init(nx, ny, "z_2");
            c_diri.init(nx, ny - 1, "c_diri");
            y_diri.init(nx, ny, "y_diri");
            for (int i = 0; i < nx; i++)
            {
                chasing_periodic_precompute(
                    x_diag[i], ny, z_1.get_ptr(i), z_2.get_ptr(i), c_diri.get_ptr(i));
            }
        }
        else
        {
            c.init(nx, ny - 1, "c");
            for (int i = 0; i < nx; i++)
            {
                chasing_standard_precompute(
                    x_diag[i], ny, c.get_ptr(i), BoundTypeYNegative, BoundTypeYPositive);
            }
        }
    }

    void chasing(field2& f, field2& p)
    {
        if constexpr (BoundTypeYNegative == PDEBoundaryType::Periodic &&
                      BoundTypeYPositive == PDEBoundaryType::Periodic)
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
                                                        BoundTypeYNegative,
                                                        BoundTypeYPositive);
                }
            }
        }
    }

private:
    double* x_diag          = nullptr;
    bool     is_no_Dirichlet = false;
    bool     has_last_vector = true;
    int      nx = 0, ny = 0;

    // Intermediate vectors
    field2 y;                        // For all boundary types
    field2 z_1, z_2, c_diri, y_diri; // For periodic boundary types
    field2 c;                        // For standard boundary types
};

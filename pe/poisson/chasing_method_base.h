#pragma once

#include "pch.h"

#include "core/base/location_boundary.h"
// #include "core/parallel/mpi/mpi_misc.h"

/**
 * @brief The class for chasing method in single direction.
 */
class ChasingMethodBase
{
protected:
    void chasing_standard_without_precompute(double            x_diag_single,
                                             int               n,
                                             double*           f,
                                             double*           c,
                                             double*           y,
                                             double*           x,
                                             PDEBoundaryType BoundaryTypeStart,
                                             PDEBoundaryType BoundaryTypeEnd);
    void chasing_standard_with_precompute(double            x_diag_single,
                                          int               n,
                                          double*           f,
                                          double*           c,
                                          double*           y,
                                          double*           x,
                                          PDEBoundaryType BoundaryTypeStart,
                                          PDEBoundaryType BoundaryTypeEnd);
    void chasing_standard_precompute(double            x_diag_single,
                                     int               n,
                                     double*           c,
                                     PDEBoundaryType BoundaryTypeStart,
                                     PDEBoundaryType BoundaryTypeEnd);
    void chasing_standard_singular(int n, double* f, double* x);
    void chasing_periodic_with_precompute(double  x_diag_single,
                                          int     n,
                                          double* f,
                                          double* c_diri,
                                          double* y_diri,
                                          double* z_1,
                                          double* z_2,
                                          double* y,
                                          double* x);
    void chasing_periodic_precompute(double x_diag_single, int n, double* z_1, double* z_2, double* c_diri);
    void chasing_periodic_singular(int n, double* f, double* x);

private:
    double p_const = 0.0;
};
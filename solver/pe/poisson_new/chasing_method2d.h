#pragma once

#include "pch.h"

#include "core/boundary/boundary_type.h"
#include "core/parallel/mpi/mpi_misc.h"

#include "chasing_method_base.h"

class ChasingMethod2D : public ChasingMethodBase
{
public:
    void init(int nx_in, int ny_in, double* _x_diag, bool _is_no_Dirichlet, bool _has_last_vector,
              PDEBoundaryType boundaryTypeYNegative, PDEBoundaryType boundaryTypeYPositive);

    void chasing(field2& f, field2& p);

private:
    double*       x_diag          = nullptr;
    bool          is_no_Dirichlet = false;
    bool          has_last_vector = true;
    int           nx = 0, ny = 0;
    PDEBoundaryType BoundaryTypeYNegative;
    PDEBoundaryType BoundaryTypeYPositive;

    // Intermediate vectors
    field2 y;                        // For all boundary types
    field2 z_1, z_2, c_diri, y_diri; // For periodic boundary types
    field2 c;                        // For standard boundary types
};

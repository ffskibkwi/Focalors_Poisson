#pragma once

#include "base/pch.h"

#include "base/location_boundary.h"
#include "base/parallel/mpi/mpi_misc.h"

#include "chasing_method_base.h"

class ChasingMethod3D : public ChasingMethodBase
{
public:
    void init(int             nx_in,
              int             ny_in,
              int             nz_in,
              double**        _x_diag,
              bool            _is_no_Dirichlet,
              bool            _has_last_vector,
              PDEBoundaryType boundary_type_down,
              PDEBoundaryType boundary_type_up);

    void chasing(field3& f, field3& p);

private:
    double**        x_diag          = nullptr;
    bool            is_no_Dirichlet = false;
    bool            has_last_vector = true;
    int             nx = 0, ny = 0, nz = 0;
    PDEBoundaryType boundary_type_down;
    PDEBoundaryType boundary_type_up;

    // Intermediate vectors
    field3 y;                        // For all boundary types
    field3 z_1, z_2, c_diri, y_diri; // For periodic boundary types
    field3 c;                        // For standard boundary types
};

#include "core/domain/domain2d.h"
#include "core/base/field2.h"
#include "core/base/location_boundary.h"
#include "pe/poisson_new/poisson_solver2d.h"

#include "io/csv_writer_2d.h"

int main(int argc, char* argv[])
{
    Domain2DUniform Omega(20, 20, 1.0, 1.0, "Omega");
    Omega.set_boundary(LocationType::Up,    PDEBoundaryType::Dirichlet);
    Omega.set_boundary(LocationType::Down,    PDEBoundaryType::Dirichlet);
    Omega.set_boundary(LocationType::Left,    PDEBoundaryType::Dirichlet);
    Omega.set_boundary(LocationType::Right,    PDEBoundaryType::Dirichlet);

    field2 Omega_field;
    Omega.construct_field(Omega_field);

    for (int i = 3; i < Omega_field.get_nx() - 3; i++)
    {
        for (int j = 3; j < Omega_field.get_ny() - 3; j++)
        {
            Omega_field(i, j) = 1.0;
        }
    }

    IO::field_to_csv(Omega_field, "result/before.txt");

    PoissonSolver2D solver(&Omega);
    solver.solve(Omega_field);

    IO::field_to_csv(Omega_field, "result/after.txt");

    field2 Omega_field_valid(Omega_field);
    for (int i = 0; i < Omega_field_valid.get_nx(); i++)
    {
        for (int j = 0; j < Omega_field_valid.get_ny(); j++)
        {
            Omega_field_valid(i, j) = -4.0 * Omega_field(i, j);
            if (i > 0)
                Omega_field_valid(i, j) += Omega_field(i - 1, j);
            if (i < Omega_field_valid.get_nx() - 1)
                Omega_field_valid(i, j) += Omega_field(i + 1, j);
            if (j > 0)
                Omega_field_valid(i, j) += Omega_field(i, j - 1);
            if (j < Omega_field_valid.get_ny() - 1)
                Omega_field_valid(i, j) += Omega_field(i, j + 1);
            Omega_field_valid(i, j) /= Omega.hx * Omega.hx;
        }
    }

    IO::field_to_csv(Omega_field_valid, "result/valid.txt");

    return 0;
}


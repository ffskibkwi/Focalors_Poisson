#include "pe/poisson/poisson_solver2d.h"

void test(PDEBoundaryType boundary_type_left,
          PDEBoundaryType boundary_type_right,
          PDEBoundaryType boundary_type_down,
          PDEBoundaryType boundary_type_up)
{
    std::cout << "Testing " << boundary_type_left << ' ' << boundary_type_right << ' ' << boundary_type_down << ' '
              << boundary_type_up << std::endl;

    int nx = 3;
    int ny = 4;
    int hx = 1.0;
    int hy = 2.0;

    PoissonSolver2D solver(
        nx, ny, hx, hy, boundary_type_left, boundary_type_right, boundary_type_down, boundary_type_up);

    field2 f(nx, ny);
    for (int i = 0; i < f.get_nx(); i++)
    {
        for (int j = 0; j < f.get_ny(); j++)
        {
            f(i, j) = i * f.get_ny() + j + 1;
        }
    }

    solver.solve(f);

    f.print();
}

// Count the number of set bits (popcount) in an integer
int countNeumanns(int mask)
{
    int count = 0;
    while (mask)
    {
        count += mask & 1;
        mask >>= 1;
    }
    return count;
}

int main()
{
    const int numParams         = 4;
    const int totalCombinations = 1 << numParams; // 2^4 = 16

    // Iterate over all possible numbers of Neumann boundaries: 4, 3, 2, 1, 0
    for (int targetNeumannCount = numParams; targetNeumannCount >= 0; --targetNeumannCount)
    {
        // Check every possible combination (represented as a bitmask)
        for (int mask = 0; mask < totalCombinations; ++mask)
        {
            if (countNeumanns(mask) == targetNeumannCount)
            {
                // Convert bitmask to boundary types:
                // bit 0 -> first parameter, bit 1 -> second, etc.
                PDEBoundaryType args[numParams];
                for (int i = 0; i < numParams; ++i)
                {
                    bool isNeumann = (mask >> i) & 1;
                    args[i]        = isNeumann ? PDEBoundaryType::Neumann : PDEBoundaryType::Dirichlet;
                }
                test(args[0], args[1], args[2], args[3]);
            }
        }
    }

    // Case 1: All four are Periodic
    test(PDEBoundaryType::Periodic, PDEBoundaryType::Periodic, PDEBoundaryType::Periodic, PDEBoundaryType::Periodic);

    // Case 2: First two are Periodic, last two vary over {Dirichlet, Neumann}
    for (int mask = 0; mask < 4; ++mask)
    {                                             // 2 bits => 4 combinations
        bool            isNeu2 = (mask >> 0) & 1; // third parameter (index 2)
        bool            isNeu3 = (mask >> 1) & 1; // fourth parameter (index 3)
        PDEBoundaryType b3     = isNeu2 ? PDEBoundaryType::Neumann : PDEBoundaryType::Dirichlet;
        PDEBoundaryType b4     = isNeu3 ? PDEBoundaryType::Neumann : PDEBoundaryType::Dirichlet;
        test(PDEBoundaryType::Periodic, PDEBoundaryType::Periodic, b3, b4);
    }

    // Case 3: Last two are Periodic, first two vary over {Dirichlet, Neumann}
    for (int mask = 0; mask < 4; ++mask)
    {
        bool            isNeu0 = (mask >> 0) & 1; // first parameter (index 0)
        bool            isNeu1 = (mask >> 1) & 1; // second parameter (index 1)
        PDEBoundaryType b1     = isNeu0 ? PDEBoundaryType::Neumann : PDEBoundaryType::Dirichlet;
        PDEBoundaryType b2     = isNeu1 ? PDEBoundaryType::Neumann : PDEBoundaryType::Dirichlet;
        test(b1, b2, PDEBoundaryType::Periodic, PDEBoundaryType::Periodic);
    }

    return 0;
}
#include "domain2d.h"

Domain2DUniform::Domain2DUniform() = default;

Domain2DUniform::Domain2DUniform(const std::string& in_name)
    : name(in_name)
{}

Domain2DUniform::Domain2DUniform(int in_nx, int in_ny, double in_lx, double in_ly, const std::string& in_name)
    : nx(in_nx)
    , ny(in_ny)
    , lx(in_lx)
    , ly(in_ly)
    , name(in_name)
{
    hx = in_lx / in_nx;
    hy = in_ly / in_ny;
}

Domain2DUniform::~Domain2DUniform() = default;

/**
 * @brief Set the number of nodes along x-direction.
 * @param in_nx Number of nodes along x (>= 1).
 */
void Domain2DUniform::set_nx(int in_nx)
{
    nx = in_nx;
    hx = lx / in_nx;
}

/**
 * @brief Set the number of nodes along y-direction.
 * @param in_ny Number of nodes along y (>= 1).
 */
void Domain2DUniform::set_ny(int in_ny)
{
    ny = in_ny;
    hy = ly / in_ny;
}

/**
 * @brief Set the domain size along x-direction.
 * @param in_lx Physical length along x (> 0).
 */
void Domain2DUniform::set_lx(double in_lx)
{
    lx = in_lx;
    hx = in_lx / nx;
}

/**
 * @brief Set the domain size along y-direction.
 * @param in_ly Physical length along y (> 0).
 */
void Domain2DUniform::set_ly(double in_ly)
{
    ly = in_ly;
    hy = in_ly / ny;
}

/**
 * @brief Set the domain size along both directions.
 * @param in_lx Physical length along x (> 0).
 * @param in_ly Physical length along y (> 0).
 */
void Domain2DUniform::set_size(double in_lx, double in_ly)
{
    lx = in_lx;
    ly = in_ly;
    hx = in_lx / nx;
    hy = in_ly / ny;
}

/**
 * @brief Set boundary condition on a given side.
 * @param loc Location side where the boundary condition applies.
 * @param type Boundary type to set.
 */
void Domain2DUniform::set_boundary(LocationType loc, PDEBoundaryType type)
{
    switch (loc)
    {
    case LocationType::Left:
        BoundaryTypeXNegative = type;
        break;
    case LocationType::Right:
        BoundaryTypeXPositive = type;
        break;
    case LocationType::Down:
        BoundaryTypeYNegative = type;
        break;
    case LocationType::Up:
        BoundaryTypeYPositive = type;
        break;
    default:
        throw std::invalid_argument("Unsupported LocationType for 2D domain");
    }
}

/**
 * @brief Get boundary condition on a given side.
 * @param loc Location side where the boundary condition applies.
 * @return Boundary type.
 */
PDEBoundaryType Domain2DUniform::get_boundary(LocationType loc) const
{
    switch (loc)
    {
    case LocationType::Left:
        return BoundaryTypeXNegative;
    case LocationType::Right:
        return BoundaryTypeXPositive;
    case LocationType::Down:
        return BoundaryTypeYNegative;
    case LocationType::Up:
        return BoundaryTypeYPositive;
    default:
        throw std::invalid_argument("Unsupported LocationType for 2D domain");
    }
}

/**
 * @brief Check whether domain profile is valid.
 * @return true if nx, ny, lx and ly are all positive; false otherwise.
 */
bool Domain2DUniform::check_profile() const { return nx > 0 && ny > 0 && lx > 0.0 && ly > 0.0; }

/**
 * @brief Check whether domain boundary is valid.
 * @return true if all boundaries are not Null; false otherwise.
 */
bool Domain2DUniform::check_boundary() const { return BoundaryTypeXNegative != PDEBoundaryType::Null && BoundaryTypeXPositive != PDEBoundaryType::Null && BoundaryTypeYNegative != PDEBoundaryType::Null && BoundaryTypeYPositive != PDEBoundaryType::Null; }



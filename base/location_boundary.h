#pragma once

#include <cstdint>
#include <stdexcept>

enum class LocationType : std::uint8_t
{
    Left,
    Right,
    Down,
    Up,
    Front,
    Back
};

/**
 * @brief Get the opposite location type
 * @param location The location type
 * @return The opposite location type
 */
inline LocationType opposite(LocationType location)
{
    switch (location)
    {
    case LocationType::Left:
        return LocationType::Right;
    case LocationType::Right:
        return LocationType::Left;
    case LocationType::Down:
        return LocationType::Up;
    case LocationType::Up:
        return LocationType::Down;
    case LocationType::Front:
        return LocationType::Back;
    case LocationType::Back:
        return LocationType::Front;
    default:
        throw std::invalid_argument("Invalid location type");
    }
}

enum class PDEBoundaryType : std::uint8_t
{
    Dirichlet,
    Neumann,
    Periodic,
    Null        //Default boundary type
};

enum class FluidBoundaryType : std::uint8_t
{
    Wall,       //Wall boundary; velocity -> Dirichlet; pressure -> Neumann
    Far,        //Far-field boundary; velocity -> Neumann; pressure -> Neumann
    Periodic,   //Periodic boundary; velocity -> Periodic; pressure -> Periodic
    Convective, //Convective boundary; velocity -> Dirichlet; pressure -> Neumann (Not sure)
    Null        //Default boundary type
};

enum class ParticleBoundaryType : std::uint8_t
{
    Bounce,
    Deposition,
    Periodic
};

enum class VariableLocationType : std::uint8_t
{
    Center,
    XEdge,
    YEdge,
    Null
};
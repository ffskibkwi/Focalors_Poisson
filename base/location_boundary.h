#pragma once

#include <cstdint>
#include <functional>
#include <iostream>
#include <stdexcept>

enum class LocationType : std::uint8_t
{
    XNegative,
    XPositive,
    YNegative,
    YPositive,
    ZNegative,
    ZPositive
};

constexpr std::array<LocationType, 4> kBoundaryLocations2D = {
    LocationType::XNegative,
    LocationType::XPositive,
    LocationType::YNegative,
    LocationType::YPositive,
};

constexpr std::array<LocationType, 6> kBoundaryLocations3D = {
    LocationType::XNegative,
    LocationType::XPositive,
    LocationType::YNegative,
    LocationType::YPositive,
    LocationType::ZNegative,
    LocationType::ZPositive,
};

std::ostream& operator<<(std::ostream& os, LocationType type);

/**
 * @brief Get the opposite location type
 * @param location The location type
 * @return The opposite location type
 */
inline LocationType opposite(LocationType location)
{
    switch (location)
    {
        case LocationType::XNegative:
            return LocationType::XPositive;
        case LocationType::XPositive:
            return LocationType::XNegative;
        case LocationType::YNegative:
            return LocationType::YPositive;
        case LocationType::YPositive:
            return LocationType::YNegative;
        case LocationType::ZNegative:
            return LocationType::ZPositive;
        case LocationType::ZPositive:
            return LocationType::ZNegative;
        default:
            throw std::invalid_argument("Invalid location type");
    }
}

enum class PDEBoundaryType : std::uint8_t
{
    Dirichlet,
    Neumann,
    Periodic,
    Adjacented,
    Null // Default boundary type
};

std::ostream& operator<<(std::ostream& os, PDEBoundaryType type);

bool isDirLike(PDEBoundaryType t);

enum class FluidBoundaryType : std::uint8_t
{
    Wall,       // Wall boundary; velocity -> Dirichlet; pressure -> Neumann
    Far,        // Far-field boundary; velocity -> Neumann; pressure -> Neumann
    Periodic,   // Periodic boundary; velocity -> Periodic; pressure -> Periodic
    Convective, // Convective boundary; velocity -> Dirichlet; pressure -> Neumann (Not sure)
    Null        // Default boundary type
};

enum class ParticleBoundaryType : std::uint8_t
{
    Bounce,
    Deposition,
    Periodic
};

enum class VariablePositionType : std::uint8_t
{
    Center,
    XFace,
    YFace,
    ZFace,
    Corner,
    Null
};

std::ostream& operator<<(std::ostream& os, VariablePositionType type);

namespace std
{
    template<>
    struct hash<LocationType>
    {
        std::size_t operator()(const LocationType& k) const
        {
            return std::hash<std::uint8_t>()(static_cast<std::uint8_t>(k));
        }
    };
} // namespace std
#pragma once

#include <cstdint>

enum class LocationType : std::uint8_t
{
    Left,
    Right,
    Down,
    Up,
    Front,
    Back
};

enum class FluidBoundaryType : std::uint8_t
{
    Dirichlet,
    Neumann,
    Periodic,
    Convective,
};

enum class ParticleBoundaryType : std::uint8_t
{
    Bounce,
    Deposition,
    Periodic,
};
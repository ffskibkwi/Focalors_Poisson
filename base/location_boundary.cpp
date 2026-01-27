#include "location_boundary.h"

std::ostream& operator<<(std::ostream& os, LocationType type)
{
    switch (type)
    {
        case LocationType::Left:
            os << "Left";
            break;
        case LocationType::Right:
            os << "Right";
            break;
        case LocationType::Down:
            os << "Down";
            break;
        case LocationType::Up:
            os << "Up";
            break;
        case LocationType::Front:
            os << "Front";
            break;
        case LocationType::Back:
            os << "Back";
            break;
        default:
            os << "Unknown";
            break;
    }
    return os;
}

std::ostream& operator<<(std::ostream& os, PDEBoundaryType type)
{
    switch (type)
    {
        case PDEBoundaryType::Dirichlet:
            os << "Dirichlet";
            break;
        case PDEBoundaryType::Neumann:
            os << "Neumann";
            break;
        case PDEBoundaryType::Periodic:
            os << "Periodic";
            break;
        case PDEBoundaryType::Adjacented:
            os << "Adjacented";
            break;
        case PDEBoundaryType::Null:
            os << "Null";
            break;
        default:
            os << "Unknown";
            break;
    }
    return os;
}

bool isDirLike(PDEBoundaryType t) { return t == PDEBoundaryType::Dirichlet || t == PDEBoundaryType::Adjacented; }
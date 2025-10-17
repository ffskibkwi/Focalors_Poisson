#pragma once

#include "base/pch.h"
#include "base/location_boundary.h"
#include "base/mesh_profile/mesh_profile_2d_base.h"
#include <stdexcept>

class Geometry2D;

class Domain2DUniform
{
public:
    // Number of nodes on each dimension
    int    nx = 0;
    int    ny = 0;
    // Size of the domain
    double lx = 0.0;
    double ly = 0.0;
    // Spatial step
    double hx = 0.0;
    double hy = 0.0;

    std::string name;
    Geometry2D* parent = nullptr;

    // Boundary conditions (four directions) and setting flags
    PDEBoundaryType BoundaryTypeXNegative = PDEBoundaryType::Null; // Left
    PDEBoundaryType BoundaryTypeXPositive = PDEBoundaryType::Null; // Right
    PDEBoundaryType BoundaryTypeYNegative = PDEBoundaryType::Null; // Down
    PDEBoundaryType BoundaryTypeYPositive = PDEBoundaryType::Null; // Up

    Domain2DUniform();
    Domain2DUniform(const std::string& in_name);
    Domain2DUniform(int in_nx, int in_ny, double in_lx, double in_ly, const std::string& in_name);

    ~Domain2DUniform();

    void set_nx(int in_nx);
    void set_ny(int in_ny);
    void set_lx(double in_lx);
    void set_ly(double in_ly);
    void set_size(double in_lx, double in_ly);

    void set_boundary(LocationType loc, PDEBoundaryType type);
    PDEBoundaryType get_boundary(LocationType loc) const;

    bool check_profile() const;
    bool check_boundary() const;

    void construct_field(field2& f);
};
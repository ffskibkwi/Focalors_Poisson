#pragma once

#include "base/field/field2.h"
#include "base/location_boundary.h"
#include "base/mesh_profile/mesh_profile_2d_base.h"

#include <stdexcept>

class Geometry2D;

class Domain2DUniform
{
public:
    // Number of nodes on each dimension
    int nx = 0;
    int ny = 0;
    // Size of the domain
    double lx = 0.0;
    double ly = 0.0;
    // Spatial step
    double hx = 0.0;
    double hy = 0.0;

    // Basic position (the position of the left-down corner)
    double pos_x, pos_y;

    std::string name;
    Geometry2D* parent = nullptr;

    // Boundary conditions (four directions) and setting flags
    PDEBoundaryType boundary_type_left  = PDEBoundaryType::Null; // Left
    PDEBoundaryType boundary_type_right = PDEBoundaryType::Null; // Right
    PDEBoundaryType boundary_type_down  = PDEBoundaryType::Null; // Down
    PDEBoundaryType boundary_type_up    = PDEBoundaryType::Null; // Up

    Domain2DUniform();
    Domain2DUniform(const std::string& in_name);
    Domain2DUniform(int in_nx, int in_ny, double in_lx, double in_ly, const std::string& in_name);
    Domain2DUniform(int in_nx, int in_ny, const std::string& in_name);

    ~Domain2DUniform();

    void set_nx(int in_nx);
    void set_ny(int in_ny);
    void set_lx(double in_lx);
    void set_ly(double in_ly);
    void set_spatial_step(double hx, double hy);
    void set_size(double in_lx, double in_ly);
    void set_position(double in_pos_x, double in_pos_y);

    double get_pos_x() const;
    double get_pos_y() const;
    double get_hx() const;
    double get_hy() const;
    double get_lx() const;
    double get_ly() const;
    int    get_nx() const;
    int    get_ny() const;

    // void set_boundary(LocationType loc, PDEBoundaryType type);
    // PDEBoundaryType get_boundary(LocationType loc) const;

    bool check_profile() const;
    bool check_boundary() const;

    void construct_field(field2& f);
};
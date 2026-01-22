#pragma once

#include "base/field/field3.h"
#include "base/location_boundary.h"

#include <stdexcept>

class Geometry3D;

class Domain3DUniform
{
public:
    // Number of nodes on each dimension
    int nx = 0;
    int ny = 0;
    int nz = 0;
    // Size of the domain
    double lx = 0.0;
    double ly = 0.0;
    double lz = 0.0;
    // Spatial step
    double hx = 0.0;
    double hy = 0.0;
    double hz = 0.0;

    // Basic position (the position of the left-down corner)
    double pos_x, pos_y, pos_z;

    // Global offset
    double offset_x = 0.0;
    double offset_y = 0.0;
    double offset_z = 0.0;

    std::string name;
    Geometry3D* parent = nullptr;

    // Boundary conditions (four directions) and setting flags
    PDEBoundaryType boundary_type_left  = PDEBoundaryType::Null; // Left
    PDEBoundaryType boundary_type_right = PDEBoundaryType::Null; // Right
    PDEBoundaryType boundary_type_front = PDEBoundaryType::Null; // Front
    PDEBoundaryType boundary_type_back  = PDEBoundaryType::Null; // Back
    PDEBoundaryType boundary_type_down  = PDEBoundaryType::Null; // Down
    PDEBoundaryType boundary_type_up    = PDEBoundaryType::Null; // Up

    Domain3DUniform();
    Domain3DUniform(const std::string& in_name);
    Domain3DUniform(int in_nx, int in_ny, int in_nz, const std::string& in_name);
    Domain3DUniform(int                in_nx,
                    int                in_ny,
                    int                in_nz,
                    double             in_lx,
                    double             in_ly,
                    double             in_lz,
                    const std::string& in_name);

    ~Domain3DUniform();

    void set_nx(int in_nx);
    void set_ny(int in_ny);
    void set_nz(int in_nz);
    void set_lx(double in_lx);
    void set_ly(double in_ly);
    void set_lz(double in_lz);
    void set_spatial_step(double hx, double hy, double hz);
    void set_size(double in_lx, double in_ly, double in_lz);
    void set_position(double in_pos_x, double in_pos_y, double in_pos_z);

    void   set_offset_x(double in_offset_x);
    void   set_offset_y(double in_offset_y);
    void   set_offset_z(double in_offset_z);
    double get_offset_x() const;
    double get_offset_y() const;
    double get_offset_z() const;

    double get_pos_x() const;
    double get_pos_y() const;
    double get_pos_z() const;
    double get_hx() const;
    double get_hy() const;
    double get_hz() const;
    double get_lx() const;
    double get_ly() const;
    double get_lz() const;
    int    get_nx() const;
    int    get_ny() const;
    int    get_nz() const;

    // void set_boundary(LocationType loc, PDEBoundaryType type);
    // PDEBoundaryType get_boundary(LocationType loc) const;

    bool check_profile() const;
    bool check_boundary() const;

    void construct_field(field3& f);
};

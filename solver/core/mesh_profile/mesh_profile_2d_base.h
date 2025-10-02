#pragma once

class MeshProfile2DBase
{
public:
    int    nx, ny;
    double lx, ly;
    bool is_uniform;

    MeshProfile2DBase(int nx, int ny, double lx, double ly)
        : nx(nx)
        , ny(ny)
        , lx(lx)
        , ly(ly)
    {}
    ~MeshProfile2DBase() = default;

    // For concat geometry, or MPI environment

    double start_x = 0.0; // The x-position of the point at southwestern corner
    double start_y = 0.0; // The y-position of the point at southeastern corner

    virtual double get_cell_pos_x(int i)    = 0;
    virtual double get_cell_pos_y(int j)    = 0;
    virtual double get_surface_pos_x(int i) = 0;
    virtual double get_surface_pos_y(int j) = 0;
};

#define EXPOSE_MESH_PROFILE_2D_BASE(mesh_profile) \
    const int    nx = mesh_profile.nx;            \
    const int    ny = mesh_profile.ny;            \
    const double lx = mesh_profile.lx;            \
    const double ly = mesh_profile.ly;

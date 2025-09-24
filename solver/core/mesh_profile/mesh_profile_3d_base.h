#pragma once

class MeshProfile3DBase
{
public:
    int    nx, ny, nz;
    double lx, ly, lz;

    MeshProfile3DBase(int nx, int ny, int nz, double lx, double ly, double lz)
        : nx(nx)
        , ny(ny)
        , nz(nz)
        , lx(lx)
        , ly(ly)
        , lz(lz)
    {}
    ~MeshProfile3DBase() = default;

    // For concat geometry, or MPI environment

    double start_x = 0.0; // The x-position of the point at southwestern corner
    double start_y = 0.0; // The y-position of the point at southeastern corner
    double start_z = 0.0; // The z-position of the point at southwestern corner

    virtual double get_cell_pos_x(int i)    = 0;
    virtual double get_cell_pos_y(int j)    = 0;
    virtual double get_cell_pos_z(int k)    = 0;
    virtual double get_surface_pos_x(int i) = 0;
    virtual double get_surface_pos_y(int j) = 0;
    virtual double get_surface_pos_z(int k) = 0;
};

#define EXPOSE_MESH_PROFILE_3D_BASE(mesh_profile) \
    const int    nx = mesh_profile.nx;            \
    const int    ny = mesh_profile.ny;            \
    const int    nz = mesh_profile.nz;            \
    const double lx = mesh_profile.lx;            \
    const double ly = mesh_profile.ly;            \
    const double lz = mesh_profile.lz;

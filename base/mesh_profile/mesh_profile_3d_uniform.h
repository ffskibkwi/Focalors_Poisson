#pragma once

#include "mesh_profile_3d_base.h"

class MeshProfile3DUniform : public MeshProfile3DBase
{
public:
    double hx;
    double hy;
    double hz;
    double inv_hx;  // 1/hx
    double inv_hy;  // 1/hy
    double inv_hz;  // 1/hz
    double inv_hx2; // 1/(hx^2)
    double inv_hy2; // 1/(hy^2)
    double inv_hz2; // 1/(hz^2)

    MeshProfile3DUniform(int nx, int ny, int nz, double lx, double ly, double lz)
        : MeshProfile3DBase(nx, ny, nz, lx, ly, lz)
        , hx(lx / nx)
        , hy(ly / ny)
        , hz(lz / nz)
        , inv_hx(1.0 / hx)
        , inv_hy(1.0 / hy)
        , inv_hz(1.0 / hz)
        , inv_hx2(1.0 / (hx * hx))
        , inv_hy2(1.0 / (hy * hy))
        , inv_hz2(1.0 / (hz * hz))
    {}
    ~MeshProfile3DUniform() = default;

    double get_cell_pos_x(int i) override { return start_x + (i + 0.5) * hx; }
    double get_cell_pos_y(int j) override { return start_y + (j + 0.5) * hy; }
    double get_cell_pos_z(int k) override { return start_z + (k + 0.5) * hz; }
    double get_surface_pos_x(int i) override { return start_x + i * hx; }
    double get_surface_pos_y(int j) override { return start_y + j * hy; }
    double get_surface_pos_z(int k) override { return start_z + k * hz; }
};

#define EXPOSE_MESH_PROFILE_3D_UNIFORM(mesh_profile) \
    EXPOSE_MESH_PROFILE_3D_BASE(mesh_profile)        \
    const double hx      = mesh_profile.hx;          \
    const double hy      = mesh_profile.hy;          \
    const double hz      = mesh_profile.hz;          \
    const double inv_hx  = mesh_profile.inv_hx;      \
    const double inv_hy  = mesh_profile.inv_hy;      \
    const double inv_hz  = mesh_profile.inv_hz;      \
    const double inv_hx2 = mesh_profile.inv_hx2;     \
    const double inv_hy2 = mesh_profile.inv_hy2;     \
    const double inv_hz2 = mesh_profile.inv_hz2;

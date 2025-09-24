#pragma once

#include "mesh_profile_2d_base.h"

class MeshProfile2DUniform : public MeshProfile2DBase
{
public:
    double hx;
    double hy;
    double inv_hx;  // 1/hx
    double inv_hy;  // 1/hy
    double inv_hx2; // 1/(hx^2)
    double inv_hy2; // 1/(hy^2)

    MeshProfile2DUniform(int nx, int ny, double lx, double ly)
        : MeshProfile2DBase(nx, ny, lx, ly)
        , hx(lx / nx)
        , hy(ly / ny)
        , inv_hx(1.0 / hx)
        , inv_hy(1.0 / hy)
        , inv_hx2(1.0 / (hx * hx))
        , inv_hy2(1.0 / (hy * hy))
    {}
    ~MeshProfile2DUniform() = default;

    double get_cell_pos_x(int i) override { return start_x + (i + 0.5) * hx; }
    double get_cell_pos_y(int j) override { return start_y + (j + 0.5) * hy; }
    double get_surface_pos_x(int i) override { return start_x + i * hx; }
    double get_surface_pos_y(int j) override { return start_y + j * hy; }
};

#define EXPOSE_MESH_PROFILE_2D_UNIFORM(mesh_profile) \
    EXPOSE_MESH_PROFILE_2D_BASE(mesh_profile)        \
    const double hx      = mesh_profile.hx;          \
    const double hy      = mesh_profile.hy;          \
    const double inv_hx  = mesh_profile.inv_hx;      \
    const double inv_hy  = mesh_profile.inv_hy;      \
    const double inv_hx2 = mesh_profile.inv_hx2;     \
    const double inv_hy2 = mesh_profile.inv_hy2;

#include "domain3d.h"

Domain3DUniform::Domain3DUniform() = default;

Domain3DUniform::Domain3DUniform(const std::string& in_name)
    : name(in_name)
{}

Domain3DUniform::Domain3DUniform(int                in_nx,
                                 int                in_ny,
                                 int                in_nz,
                                 double             in_lx,
                                 double             in_ly,
                                 double             in_lz,
                                 const std::string& in_name)
    : nx(in_nx)
    , ny(in_ny)
    , nz(in_nz)
    , lx(in_lx)
    , ly(in_ly)
    , lz(in_lz)
    , name(in_name)
{
    hx = in_lx / in_nx;
    hy = in_ly / in_ny;
    hz = in_lz / in_nz;
}

Domain3DUniform::~Domain3DUniform() = default;

/**
 * @brief Set the number of nodes along x-direction.
 * @param in_nx Number of nodes along x (>= 1).
 */
void Domain3DUniform::set_nx(int in_nx)
{
    nx = in_nx;
    hx = lx / in_nx;
}

/**
 * @brief Set the number of nodes along y-direction.
 * @param in_ny Number of nodes along y (>= 1).
 */
void Domain3DUniform::set_ny(int in_ny)
{
    ny = in_ny;
    hy = ly / in_ny;
}

/**
 * @brief Set the number of nodes along z-direction.
 * @param in_nz Number of nodes along z (>= 1).
 */
void Domain3DUniform::set_nz(int in_nz)
{
    nz = in_nz;
    hz = lz / in_nz;
}

/**
 * @brief Set the domain size along x-direction.
 * @param in_lx Physical length along x (> 0).
 */
void Domain3DUniform::set_lx(double in_lx)
{
    lx = in_lx;
    hx = in_lx / nx;
}

/**
 * @brief Set the domain size along y-direction.
 * @param in_ly Physical length along y (> 0).
 */
void Domain3DUniform::set_ly(double in_ly)
{
    ly = in_ly;
    hy = in_ly / ny;
}

/**
 * @brief Set the domain size along z-direction.
 * @param in_lz Physical length along z (> 0).
 */
void Domain3DUniform::set_lz(double in_lz)
{
    lz = in_lz;
    hz = in_lz / nz;
}

/**
 * @brief Set the domain size along both directions.
 * @param in_lx Physical length along x (> 0).
 * @param in_ly Physical length along y (> 0).
 * @param in_lz Physical length along z (> 0).
 */
void Domain3DUniform::set_size(double in_lx, double in_ly, double in_lz)
{
    lx = in_lx;
    ly = in_ly;
    lz = in_lz;
    if (nx != 0)
        hx = in_lx / nx;
    if (ny != 0)
        hy = in_ly / ny;
    if (nz != 0)
        hz = in_lz / nz;
}

void Domain3DUniform::set_position(double in_pos_x, double in_pos_y, double in_pos_z)
{
    pos_x = in_pos_x;
    pos_y = in_pos_y;
    pos_z = in_pos_z;
}

/**
 * @brief Check whether domain profile is valid.
 * @return true if nx, ny, nz, lx, ly, lz are all positive; false otherwise.
 */
bool Domain3DUniform::check_profile() const { return nx > 0 && ny > 0 && nz > 0 && lx > 0.0 && ly > 0.0 && lz > 0.0; }

/**
 * @brief Check whether domain boundary is valid.
 * @return true if all boundaries are not Null; false otherwise.
 */
bool Domain3DUniform::check_boundary() const
{
    return boundary_type_left != PDEBoundaryType::Null && boundary_type_right != PDEBoundaryType::Null &&
           boundary_type_front != PDEBoundaryType::Null && boundary_type_back != PDEBoundaryType::Null &&
           boundary_type_down != PDEBoundaryType::Null && boundary_type_up != PDEBoundaryType::Null;
}

void Domain3DUniform::construct_field(field3& f) { f.init(nx, ny, nz); }

double Domain3DUniform::get_pos_x() const { return pos_x; }
double Domain3DUniform::get_pos_y() const { return pos_y; }
double Domain3DUniform::get_pos_z() const { return pos_z; }
double Domain3DUniform::get_hx() const { return hx; }
double Domain3DUniform::get_hy() const { return hy; }
double Domain3DUniform::get_hz() const { return hz; }
double Domain3DUniform::get_lx() const { return lx; }
double Domain3DUniform::get_ly() const { return ly; }
double Domain3DUniform::get_lz() const { return lz; }
int    Domain3DUniform::get_nx() const { return nx; }
int    Domain3DUniform::get_ny() const { return ny; }
int    Domain3DUniform::get_nz() const { return nz; }
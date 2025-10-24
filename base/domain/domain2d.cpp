#include "domain2d.h"

Domain2DUniform::Domain2DUniform() = default;

Domain2DUniform::Domain2DUniform(const std::string& in_name)
    : name(in_name)
{}

Domain2DUniform::Domain2DUniform(int in_nx, int in_ny, double in_lx, double in_ly, const std::string& in_name)
    : nx(in_nx)
    , ny(in_ny)
    , lx(in_lx)
    , ly(in_ly)
    , name(in_name)
{
    hx = in_lx / in_nx;
    hy = in_ly / in_ny;
}

Domain2DUniform::~Domain2DUniform() = default;

/**
 * @brief Set the number of nodes along x-direction.
 * @param in_nx Number of nodes along x (>= 1).
 */
void Domain2DUniform::set_nx(int in_nx)
{
    nx = in_nx;
    hx = lx / in_nx;
}

/**
 * @brief Set the number of nodes along y-direction.
 * @param in_ny Number of nodes along y (>= 1).
 */
void Domain2DUniform::set_ny(int in_ny)
{
    ny = in_ny;
    hy = ly / in_ny;
}

/**
 * @brief Set the domain size along x-direction.
 * @param in_lx Physical length along x (> 0).
 */
void Domain2DUniform::set_lx(double in_lx)
{
    lx = in_lx;
    hx = in_lx / nx;
}

/**
 * @brief Set the domain size along y-direction.
 * @param in_ly Physical length along y (> 0).
 */
void Domain2DUniform::set_ly(double in_ly)
{
    ly = in_ly;
    hy = in_ly / ny;
}

/**
 * @brief Set the domain size along both directions.
 * @param in_lx Physical length along x (> 0).
 * @param in_ly Physical length along y (> 0).
 */
void Domain2DUniform::set_size(double in_lx, double in_ly)
{
    lx = in_lx;
    ly = in_ly;
    if (nx != 0 )
        hx = in_lx / nx;
    if (ny != 0)
        hy = in_ly / ny;
}

void Domain2DUniform::set_position(double in_pos_x, double in_pos_y)
{
    pos_x = in_pos_x;
    pos_y = in_pos_y;
}

/**
 * @brief Check whether domain profile is valid.
 * @return true if nx, ny, lx and ly are all positive; false otherwise.
 */
bool Domain2DUniform::check_profile() const { return nx > 0 && ny > 0 && lx > 0.0 && ly > 0.0; }

/**
 * @brief Check whether domain boundary is valid.
 * @return true if all boundaries are not Null; false otherwise.
 */
bool Domain2DUniform::check_boundary() const { return boundary_type_left != PDEBoundaryType::Null && boundary_type_right != PDEBoundaryType::Null && boundary_type_down != PDEBoundaryType::Null && boundary_type_up != PDEBoundaryType::Null; }

void Domain2DUniform::construct_field(field2& f)
{
    f.init(nx, ny);
}

double Domain2DUniform::get_pos_x() const { return pos_x; }
double Domain2DUniform::get_pos_y() const { return pos_y; }
double Domain2DUniform::get_hx() const { return hx; }
double Domain2DUniform::get_hy() const { return hy; }
double Domain2DUniform::get_lx() const { return lx; }
double Domain2DUniform::get_ly() const { return ly; }
int Domain2DUniform::get_nx() const { return nx; }
int Domain2DUniform::get_ny() const { return ny; }
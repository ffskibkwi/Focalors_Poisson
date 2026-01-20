#include "variable3d.h"
#include <algorithm>

Variable3D::Variable3D(const std::string& in_name)
    : name(in_name)
{}

Variable3D::~Variable3D() { cleanup_buffers(); }

void Variable3D::cleanup_buffers()
{
    // Clean up boundary value buffers
    for (auto& domainPair : boundary_value_map)
    {
        for (auto& locPair : domainPair.second)
        {
            delete locPair.second;
        }
    }

    // Clean up buffer_map
    for (auto& domainPair : buffer_map)
    {
        for (auto& locPair : domainPair.second)
        {
            delete locPair.second;
        }
    }

    if (position_type == VariablePositionType::XFaceCenter)
    {
        for (auto& pair : corner_value_map_y)
            delete[] pair.second;
        for (auto& pair : corner_value_map_z)
            delete[] pair.second;
    }
    else if (position_type == VariablePositionType::YFaceCenter)
    {
        for (auto& pair : corner_value_map_x)
            delete[] pair.second;
        for (auto& pair : corner_value_map_z)
            delete[] pair.second;
    }
    else if (position_type == VariablePositionType::ZFaceCenter)
    {
        for (auto& pair : corner_value_map_x)
            delete[] pair.second;
        for (auto& pair : corner_value_map_y)
            delete[] pair.second;
    }
}

/**
 * @brief Attach this variable to a geometry.
 * @param g Geometry to bind to this variable.
 */
void Variable3D::set_geometry(Geometry3D& g)
{
    geometry = &g;
    for (auto& domainAdjPair : geometry->adjacency)
    {
        Domain3DUniform* domain = domainAdjPair.first;
        for (auto& locToNeighbor : domainAdjPair.second)
        {
            LocationType loc               = locToNeighbor.first;
            boundary_type_map[domain][loc] = PDEBoundaryType::Adjacented;
        }
    }
}

/**
 * @brief Check geometry (before set fields)
 * @throws std::runtime_error If geometry is not set, domain is not part of the geometry,
 *         or the domain mesh size is invalid.
 */
void Variable3D::check_geometry(Domain3DUniform* s)
{
    if (geometry == nullptr)
        throw std::runtime_error("Variable3D has no geometry set");
    if (s->parent != geometry)
        throw std::runtime_error("Domain not found in geometry");
    if (s->nx <= 0 || s->ny <= 0 || s->nz <= 0)
        throw std::runtime_error("Domain mesh size invalid");
}

/**
 * @brief Initialize a field and bind it to a domain within the attached geometry.
 * @param s Domain to which the field belongs.
 * @param f Field storage to initialize and register.
 * @throws std::runtime_error If geometry is not set, domain is not part of the geometry,
 *         or the domain mesh size is invalid.
 *
 * If the field name is "Default", it will be initialized with the name "variable".
 */
void Variable3D::set_center_field(Domain3DUniform* s, field3& f)
{
    check_geometry(s);

    if (f.get_name() == "Default")
        f.init(s->nx, s->ny, s->nz, name + "_" + s->name);
    else
        f.init(s->nx, s->ny, s->nz);

    field_map[s] = &f;

    // Center variable in 3D needs buffers for Left, Right, Front, Back, Down, Up
    // Left/Right: y-z plane (ny * nz)
    buffer_map[s][LocationType::Left]  = new field2(s->ny, s->nz);
    buffer_map[s][LocationType::Right] = new field2(s->ny, s->nz);
    // Front/Back: x-z plane (nx * nz)
    buffer_map[s][LocationType::Front] = new field2(s->nx, s->nz);
    buffer_map[s][LocationType::Back]  = new field2(s->nx, s->nz);
    // Down/Up: x-y plane (nx * ny)
    buffer_map[s][LocationType::Down] = new field2(s->nx, s->ny);
    buffer_map[s][LocationType::Up]   = new field2(s->nx, s->ny);

    position_type = VariablePositionType::Center;
}

void Variable3D::set_x_face_center_field(Domain3DUniform* s, field3& f)
{
    check_geometry(s);

    if (f.get_name() == "Default")
        f.init(s->nx, s->ny, s->nz, name + "_" + s->name);
    else
        f.init(s->nx, s->ny, s->nz);

    field_map[s] = &f;

    // Center variable in 3D needs buffers for Left, Right, Front, Back, Down, Up
    // Left/Right: y-z plane (ny * nz)
    buffer_map[s][LocationType::Left]  = new field2(s->ny, s->nz);
    buffer_map[s][LocationType::Right] = new field2(s->ny, s->nz);
    // Front/Back: x-z plane (nx * nz)
    buffer_map[s][LocationType::Front] = new field2(s->nx, s->nz);
    buffer_map[s][LocationType::Back]  = new field2(s->nx, s->nz);
    // Down/Up: x-y plane (nx * ny)
    buffer_map[s][LocationType::Down] = new field2(s->nx, s->ny);
    buffer_map[s][LocationType::Up]   = new field2(s->nx, s->ny);

    corner_value_map_y[s] = new double[s->ny];
    corner_value_map_z[s] = new double[s->nz];

    position_type = VariablePositionType::XFaceCenter;
}

void Variable3D::set_y_face_center_field(Domain3DUniform* s, field3& f)
{
    check_geometry(s);

    if (f.get_name() == "Default")
        f.init(s->nx, s->ny, s->nz, name + "_" + s->name);
    else
        f.init(s->nx, s->ny, s->nz);

    field_map[s] = &f;

    // Center variable in 3D needs buffers for Left, Right, Front, Back, Down, Up
    // Left/Right: y-z plane (ny * nz)
    buffer_map[s][LocationType::Left]  = new field2(s->ny, s->nz);
    buffer_map[s][LocationType::Right] = new field2(s->ny, s->nz);
    // Front/Back: x-z plane (nx * nz)
    buffer_map[s][LocationType::Front] = new field2(s->nx, s->nz);
    buffer_map[s][LocationType::Back]  = new field2(s->nx, s->nz);
    // Down/Up: x-y plane (nx * ny)
    buffer_map[s][LocationType::Down] = new field2(s->nx, s->ny);
    buffer_map[s][LocationType::Up]   = new field2(s->nx, s->ny);

    corner_value_map_x[s] = new double[s->nx];
    corner_value_map_z[s] = new double[s->nz];

    position_type = VariablePositionType::YFaceCenter;
}

void Variable3D::set_z_face_center_field(Domain3DUniform* s, field3& f)
{
    check_geometry(s);

    if (f.get_name() == "Default")
        f.init(s->nx, s->ny, s->nz, name + "_" + s->name);
    else
        f.init(s->nx, s->ny, s->nz);

    field_map[s] = &f;

    // Center variable in 3D needs buffers for Left, Right, Front, Back, Down, Up
    // Left/Right: y-z plane (ny * nz)
    buffer_map[s][LocationType::Left]  = new field2(s->ny, s->nz);
    buffer_map[s][LocationType::Right] = new field2(s->ny, s->nz);
    // Front/Back: x-z plane (nx * nz)
    buffer_map[s][LocationType::Front] = new field2(s->nx, s->nz);
    buffer_map[s][LocationType::Back]  = new field2(s->nx, s->nz);
    // Down/Up: x-y plane (nx * ny)
    buffer_map[s][LocationType::Down] = new field2(s->nx, s->ny);
    buffer_map[s][LocationType::Up]   = new field2(s->nx, s->ny);

    corner_value_map_x[s] = new double[s->nx];
    corner_value_map_y[s] = new double[s->ny];

    position_type = VariablePositionType::ZFaceCenter;
}

void Variable3D::set_boundary_type(Domain3DUniform* s, LocationType loc, PDEBoundaryType type)
{
    check_geometry(s);
    // 若该侧为几何邻接面，则仅允许设置为 Adjacented
    if (geometry && geometry->adjacency.count(s) && geometry->adjacency[s].count(loc))
    {
        if (type != PDEBoundaryType::Adjacented)
            throw std::runtime_error("Attempt to override an adjacented face with non-Adjacented boundary on domain " +
                                     s->name);
        boundary_type_map[s][loc] = PDEBoundaryType::Adjacented;
        return;
    }

    // 若该侧已被设置为 Adjacented，则不允许再改为其他类型
    if (boundary_type_map[s].count(loc) && boundary_type_map[s][loc] == PDEBoundaryType::Adjacented &&
        type != PDEBoundaryType::Adjacented)
    {
        throw std::runtime_error("Attempt to change previously Adjacented boundary to another type on domain " +
                                 s->name);
    }

    boundary_type_map[s][loc] = type;
}

void Variable3D::set_boundary_value(Domain3DUniform* s, LocationType loc, double in_value)
{
    check_geometry(s);
    has_boundary_value_map[s][loc] = true;

    if (loc == LocationType::Left || loc == LocationType::Right)
    {
        // Left/Right: y-z plane (ny * nz)
        int size                   = s->ny * s->nz;
        boundary_value_map[s][loc] = new field2(s->ny, s->nz);
        boundary_value_map[s][loc]->clear(in_value);
    }
    else if (loc == LocationType::Front || loc == LocationType::Back)
    {
        // Front/Back: x-z plane (nx * nz)
        boundary_value_map[s][loc] = new field2(s->nx, s->nz);
        boundary_value_map[s][loc]->clear(in_value);
    }
    else if (loc == LocationType::Down || loc == LocationType::Up)
    {
        // Down/Up: x-y plane (nx * ny)
        boundary_value_map[s][loc] = new field2(s->nx, s->ny);
        boundary_value_map[s][loc]->clear(in_value);
    }
}

void Variable3D::set_boundary_value(Domain3DUniform*                              s,
                                    LocationType                                  loc,
                                    std::function<double(double, double, double)> f)
{
    check_geometry(s);
    // TODO: Implement function-based boundary value setting for 3D
    // This would require spatial coordinates (x, y, z) for each boundary point
}

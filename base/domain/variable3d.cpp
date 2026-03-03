#include "variable3d.h"
#include "base/parallel/omp/enable_openmp.h"

#include <algorithm>
#include <array>

Variable3D::Variable3D(const std::string& in_name)
    : name(in_name)
{}

Variable3D::~Variable3D()
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

    if (position_type == VariablePositionType::XFace)
    {
        for (auto& pair : corner_y_map)
            delete[] pair.second;
        for (auto& pair : corner_z_map)
            delete[] pair.second;
    }
    else if (position_type == VariablePositionType::YFace)
    {
        for (auto& pair : corner_x_map)
            delete[] pair.second;
        for (auto& pair : corner_z_map)
            delete[] pair.second;
    }
    else if (position_type == VariablePositionType::ZFace)
    {
        for (auto& pair : corner_x_map)
            delete[] pair.second;
        for (auto& pair : corner_y_map)
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

    f.init(s->nx, s->ny, s->nz, name + "_" + s->name);

    field_map[s] = &f;

    buffer_map[s][LocationType::XNegative] = new field2(s->ny, s->nz);
    buffer_map[s][LocationType::XPositive] = new field2(s->ny, s->nz);
    buffer_map[s][LocationType::YNegative] = new field2(s->nx, s->nz);
    buffer_map[s][LocationType::YPositive] = new field2(s->nx, s->nz);
    buffer_map[s][LocationType::ZNegative] = new field2(s->nx, s->ny);
    buffer_map[s][LocationType::ZPositive] = new field2(s->nx, s->ny);

    position_type = VariablePositionType::Center;
}

void Variable3D::set_x_face_center_field(Domain3DUniform* s, field3& f)
{
    check_geometry(s);

    f.init(s->nx, s->ny, s->nz, name + "_" + s->name);

    field_map[s] = &f;

    // Center variable in 3D needs buffers for XNegative, XPositive, YNegative, YPositive, ZNegative, ZPositive
    // XNegative/XPositive: y-z plane (ny * nz)
    buffer_map[s][LocationType::XNegative] = new field2(s->ny, s->nz);
    buffer_map[s][LocationType::XPositive] = new field2(s->ny, s->nz);
    // YNegative/YPositive: x-z plane (nx * nz)
    buffer_map[s][LocationType::YNegative] = new field2(s->nx, s->nz);
    buffer_map[s][LocationType::YPositive] = new field2(s->nx, s->nz);
    // ZNegative/ZPositive: x-y plane (nx * ny)
    buffer_map[s][LocationType::ZNegative] = new field2(s->nx, s->ny);
    buffer_map[s][LocationType::ZPositive] = new field2(s->nx, s->ny);

    corner_y_map[s] = new double[s->ny];
    corner_z_map[s] = new double[s->nz];

    position_type = VariablePositionType::XFace;
}

void Variable3D::set_y_face_center_field(Domain3DUniform* s, field3& f)
{
    check_geometry(s);

    f.init(s->nx, s->ny, s->nz, name + "_" + s->name);

    field_map[s] = &f;

    // Center variable in 3D needs buffers for XNegative, XPositive, YNegative, YPositive, ZNegative, ZPositive
    // XNegative/XPositive: y-z plane (ny * nz)
    buffer_map[s][LocationType::XNegative] = new field2(s->ny, s->nz);
    buffer_map[s][LocationType::XPositive] = new field2(s->ny, s->nz);
    // YNegative/YPositive: x-z plane (nx * nz)
    buffer_map[s][LocationType::YNegative] = new field2(s->nx, s->nz);
    buffer_map[s][LocationType::YPositive] = new field2(s->nx, s->nz);
    // ZNegative/ZPositive: x-y plane (nx * ny)
    buffer_map[s][LocationType::ZNegative] = new field2(s->nx, s->ny);
    buffer_map[s][LocationType::ZPositive] = new field2(s->nx, s->ny);

    corner_x_map[s] = new double[s->nx];
    corner_z_map[s] = new double[s->nz];

    position_type = VariablePositionType::YFace;
}

void Variable3D::set_z_face_center_field(Domain3DUniform* s, field3& f)
{
    check_geometry(s);

    f.init(s->nx, s->ny, s->nz, name + "_" + s->name);

    field_map[s] = &f;

    // Center variable in 3D needs buffers for XNegative, XPositive, YNegative, YPositive, ZNegative, ZPositive
    // XNegative/XPositive: y-z plane (ny * nz)
    buffer_map[s][LocationType::XNegative] = new field2(s->ny, s->nz);
    buffer_map[s][LocationType::XPositive] = new field2(s->ny, s->nz);
    // YNegative/YPositive: x-z plane (nx * nz)
    buffer_map[s][LocationType::YNegative] = new field2(s->nx, s->nz);
    buffer_map[s][LocationType::YPositive] = new field2(s->nx, s->nz);
    // ZNegative/ZPositive: x-y plane (nx * ny)
    buffer_map[s][LocationType::ZNegative] = new field2(s->nx, s->ny);
    buffer_map[s][LocationType::ZPositive] = new field2(s->nx, s->ny);

    corner_x_map[s] = new double[s->nx];
    corner_y_map[s] = new double[s->ny];

    position_type = VariablePositionType::ZFace;
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

void Variable3D::set_boundary_type(Domain3DUniform*                                                s,
                                   std::initializer_list<std::pair<LocationType, PDEBoundaryType>> list)
{
    for (const auto& pair : list)
    {
        set_boundary_type(s, pair.first, pair.second);
    }
}

void Variable3D::set_boundary_type(PDEBoundaryType type)
{
    if (geometry == nullptr)
        throw std::runtime_error("Variable2D has no geometry set");

    for (auto* domain : geometry->domains)
    {
        auto& type_map = boundary_type_map[domain];
        for (LocationType loc : kBoundaryLocations3D)
        {
            if (type_map.find(loc) != type_map.end())
                continue;
            set_boundary_type(domain, loc, type);
        }
    }
}

void Variable3D::set_boundary_value(Domain3DUniform* s, LocationType loc, double in_value)
{
    check_geometry(s);
    has_boundary_value_map[s][loc] = true;

    if (loc == LocationType::XNegative || loc == LocationType::XPositive)
    {
        // XNegative/XPositive: y-z plane (ny * nz)
        int size                   = s->ny * s->nz;
        boundary_value_map[s][loc] = new field2(s->ny, s->nz);
        boundary_value_map[s][loc]->clear(in_value);
    }
    else if (loc == LocationType::YNegative || loc == LocationType::YPositive)
    {
        // YNegative/YPositive: x-z plane (nx * nz)
        boundary_value_map[s][loc] = new field2(s->nx, s->nz);
        boundary_value_map[s][loc]->clear(in_value);
    }
    else if (loc == LocationType::ZNegative || loc == LocationType::ZPositive)
    {
        // ZNegative/ZPositive: x-y plane (nx * ny)
        boundary_value_map[s][loc] = new field2(s->nx, s->ny);
        boundary_value_map[s][loc]->clear(in_value);
    }
}

void Variable3D::set_boundary_value(Domain3DUniform*                              s,
                                    LocationType                                  loc,
                                    std::function<double(double, double, double)> f)
{
    check_geometry(s);
    has_boundary_value_map[s][loc] = true;

    double shift_x = 0.5;
    double shift_y = 0.5;
    double shift_z = 0.5;

    switch (position_type)
    {
        case VariablePositionType::Center:
            shift_x = 0.5;
            shift_y = 0.5;
            shift_z = 0.5;
            break;
        case VariablePositionType::XFace:
            shift_x = 0.0;
            shift_y = 0.5;
            shift_z = 0.5;
            break;
        case VariablePositionType::YFace:
            shift_x = 0.5;
            shift_y = 0.0;
            shift_z = 0.5;
            break;
        case VariablePositionType::ZFace:
            shift_x = 0.5;
            shift_y = 0.5;
            shift_z = 0.0;
            break;
        case VariablePositionType::Corner:
            shift_x = 0.0;
            shift_y = 0.0;
            shift_z = 0.0;
            break;
        default:
            break;
    }

    switch (loc)
    {
        case LocationType::XNegative:
        case LocationType::XPositive:
            shift_x = 0.0;
            break;
        case LocationType::YNegative:
        case LocationType::YPositive:
            shift_y = 0.0;
            break;
        case LocationType::ZNegative:
        case LocationType::ZPositive:
            shift_z = 0.0;
            break;
        default:
            break;
    }

    if (loc == LocationType::XNegative || loc == LocationType::XPositive)
    {
        int j_size, k_size;
        if (position_type == VariablePositionType::Corner)
        {
            j_size = s->ny + 1;
            k_size = s->nz + 1;
        }
        else
        {
            j_size = s->ny;
            k_size = s->nz;
        }

        boundary_value_map[s][loc] = new field2(j_size, k_size);

        int i_idx = (loc == LocationType::XNegative) ? 0 : s->nx;

        for (int j = 0; j < j_size; j++)
        {
            for (int k = 0; k < k_size; k++)
            {
                double gx                           = s->get_offset_x() + (shift_x + i_idx) * s->get_hx();
                double gy                           = s->get_offset_y() + (shift_y + j) * s->get_hy();
                double gz                           = s->get_offset_z() + (shift_z + k) * s->get_hz();
                (*boundary_value_map[s][loc])(j, k) = f(gx, gy, gz);
            }
        }
    }
    else if (loc == LocationType::YNegative || loc == LocationType::YPositive)
    {
        int i_size, k_size;
        if (position_type == VariablePositionType::Corner)
        {
            i_size = s->nx + 1;
            k_size = s->nz + 1;
        }
        else
        {
            i_size = s->nx;
            k_size = s->nz;
        }

        boundary_value_map[s][loc] = new field2(i_size, k_size);

        int j_idx = (loc == LocationType::YNegative) ? 0 : s->ny;

        for (int i = 0; i < i_size; i++)
        {
            for (int k = 0; k < k_size; k++)
            {
                double gx                           = s->get_offset_x() + (shift_x + i) * s->get_hx();
                double gy                           = s->get_offset_y() + (shift_y + j_idx) * s->get_hy();
                double gz                           = s->get_offset_z() + (shift_z + k) * s->get_hz();
                (*boundary_value_map[s][loc])(i, k) = f(gx, gy, gz);
            }
        }
    }
    else if (loc == LocationType::ZNegative || loc == LocationType::ZPositive)
    {
        int i_size, j_size;
        if (position_type == VariablePositionType::Corner)
        {
            i_size = s->nx + 1;
            j_size = s->ny + 1;
        }
        else
        {
            i_size = s->nx;
            j_size = s->ny;
        }

        boundary_value_map[s][loc] = new field2(i_size, j_size);

        int k_idx = (loc == LocationType::ZNegative) ? 0 : s->nz;

        for (int i = 0; i < i_size; i++)
        {
            for (int j = 0; j < j_size; j++)
            {
                double gx                           = s->get_offset_x() + (shift_x + i) * s->get_hx();
                double gy                           = s->get_offset_y() + (shift_y + j) * s->get_hy();
                double gz                           = s->get_offset_z() + (shift_z + k_idx) * s->get_hz();
                (*boundary_value_map[s][loc])(i, j) = f(gx, gy, gz);
            }
        }
    }
}

void Variable3D::set_boundary(std::function<double(double, double, double)> f)
{
    if (geometry == nullptr)
        throw std::runtime_error("Variable2D has no geometry set");

    for (auto kv : field_map)
    {
        Domain3DUniform* domain = kv.first;

        auto& type_map = boundary_type_map[domain];
        auto& has_map  = has_boundary_value_map[domain];
        for (LocationType loc : kBoundaryLocations3D)
        {
            auto type_it = type_map.find(loc);
            if (type_it == type_map.end())
                continue;
            if (type_it->second == PDEBoundaryType::Adjacented || type_it->second == PDEBoundaryType::Null)
                continue;
            auto has_it = has_map.find(loc);
            if (has_it != has_map.end() && has_it->second)
                continue;

            set_boundary_value(domain, loc, f);
        }
    }
}

void Variable3D::set_buffer_value(Domain3DUniform* s, LocationType loc, std::function<double(double, double, double)> f)
{
    check_geometry(s);

    double shift_x = 0.5;
    double shift_y = 0.5;
    double shift_z = 0.5;

    switch (position_type)
    {
        case VariablePositionType::Center:
            shift_x = 0.5;
            shift_y = 0.5;
            shift_z = 0.5;
            break;
        case VariablePositionType::XFace:
            shift_x = 0.0;
            shift_y = 0.5;
            shift_z = 0.5;
            break;
        case VariablePositionType::YFace:
            shift_x = 0.5;
            shift_y = 0.0;
            shift_z = 0.5;
            break;
        case VariablePositionType::ZFace:
            shift_x = 0.5;
            shift_y = 0.5;
            shift_z = 0.0;
            break;
        case VariablePositionType::Corner:
            shift_x = 0.0;
            shift_y = 0.0;
            shift_z = 0.0;
            break;
        default:
            break;
    }

    if (loc == LocationType::XNegative || loc == LocationType::XPositive)
    {
        int j_size, k_size;
        if (position_type == VariablePositionType::Corner)
        {
            j_size = s->ny + 1;
            k_size = s->nz + 1;
        }
        else
        {
            j_size = s->ny;
            k_size = s->nz;
        }

        int i_idx = (loc == LocationType::XNegative) ? -1 : s->nx;

        for (int j = 0; j < j_size; j++)
        {
            for (int k = 0; k < k_size; k++)
            {
                double gx                   = s->get_offset_x() + (shift_x + i_idx) * s->get_hx();
                double gy                   = s->get_offset_y() + (shift_y + j) * s->get_hy();
                double gz                   = s->get_offset_z() + (shift_z + k) * s->get_hz();
                (*buffer_map[s][loc])(j, k) = f(gx, gy, gz);
            }
        }
    }
    else if (loc == LocationType::YNegative || loc == LocationType::YPositive)
    {
        int i_size, k_size;
        if (position_type == VariablePositionType::Corner)
        {
            i_size = s->nx + 1;
            k_size = s->nz + 1;
        }
        else
        {
            i_size = s->nx;
            k_size = s->nz;
        }

        int j_idx = (loc == LocationType::YNegative) ? -1 : s->ny;

        for (int i = 0; i < i_size; i++)
        {
            for (int k = 0; k < k_size; k++)
            {
                double gx                   = s->get_offset_x() + (shift_x + i) * s->get_hx();
                double gy                   = s->get_offset_y() + (shift_y + j_idx) * s->get_hy();
                double gz                   = s->get_offset_z() + (shift_z + k) * s->get_hz();
                (*buffer_map[s][loc])(i, k) = f(gx, gy, gz);
            }
        }
    }
    else if (loc == LocationType::ZNegative || loc == LocationType::ZPositive)
    {
        int i_size, j_size;
        if (position_type == VariablePositionType::Corner)
        {
            i_size = s->nx + 1;
            j_size = s->ny + 1;
        }
        else
        {
            i_size = s->nx;
            j_size = s->ny;
        }

        int k_idx = (loc == LocationType::ZNegative) ? -1 : s->nz;

        for (int i = 0; i < i_size; i++)
        {
            for (int j = 0; j < j_size; j++)
            {
                double gx                   = s->get_offset_x() + (shift_x + i) * s->get_hx();
                double gy                   = s->get_offset_y() + (shift_y + j) * s->get_hy();
                double gz                   = s->get_offset_z() + (shift_z + k_idx) * s->get_hz();
                (*buffer_map[s][loc])(i, j) = f(gx, gy, gz);
            }
        }
    }
}

void Variable3D::set_corner(Domain3DUniform* s, std::function<double(double, double, double)> f)
{
    check_geometry(s);

    if (position_type == VariablePositionType::XFace)
    {
        {
            double x = s->get_offset_x() + s->nx * s->hx;
            double z = s->get_offset_z() - 0.5 * s->hz;
            for (int j = 0; j < s->ny; j++)
                corner_y_map[s][j] = f(x, (j + 0.5) * s->hy, z);
        }
        {
            double x = s->get_offset_x() + s->nx * s->hx;
            double y = s->get_offset_y() - 0.5 * s->hy;
            for (int k = 0; k < s->nz; k++)
                corner_z_map[s][k] = f(x, y, (k + 0.5) * s->hz);
        }
    }
    else if (position_type == VariablePositionType::YFace)
    {
        {
            double y = s->get_offset_y() + s->ny * s->hy;
            double z = s->get_offset_z() - 0.5 * s->hz;
            for (int i = 0; i < s->nx; i++)
                corner_x_map[s][i] = f((i + 0.5) * s->hx, y, z);
        }
        {
            double x = s->get_offset_x() - 0.5 * s->hx;
            double y = s->get_offset_y() + s->ny * s->hy;
            for (int k = 0; k < s->nz; k++)
                corner_z_map[s][k] = f(x, y, (k + 0.5) * s->hz);
        }
    }
    else if (position_type == VariablePositionType::ZFace)
    {
        {
            double y = s->get_offset_y() - 0.5 * s->hy;
            double z = s->get_offset_z() + s->nz * s->hz;
            for (int i = 0; i < s->nx; i++)
                corner_x_map[s][i] = f((i + 0.5) * s->hx, y, z);
        }
        {
            double x = s->get_offset_x() - 0.5 * s->hx;
            double z = s->get_offset_z() + s->nz * s->hz;
            for (int j = 0; j < s->ny; j++)
                corner_y_map[s][j] = f(x, (j + 0.5) * s->hy, z);
        }
    }
}

void Variable3D::set_buffer(std::function<double(double, double, double)> f)
{
    if (geometry == nullptr)
        throw std::runtime_error("Variable2D has no geometry set");

    for (auto kv : field_map)
    {
        Domain3DUniform* domain = kv.first;

        auto& type_map = boundary_type_map[domain];
        for (LocationType loc : kBoundaryLocations3D)
        {
            auto type_it = type_map.find(loc);
            if (type_it == type_map.end())
                continue;
            if (type_it->second == PDEBoundaryType::Adjacented || type_it->second == PDEBoundaryType::Null)
                continue;

            set_buffer_value(domain, loc, f);
        }
    }
}

void Variable3D::set_corner(std::function<double(double, double, double)> f)
{
    for (auto kv : field_map)
        set_corner(kv.first, f);
}

void Variable3D::set_value(std::function<double(double, double, double)> func)
{
    double shift_x = 0.5;
    double shift_y = 0.5;
    double shift_z = 0.5;

    switch (position_type)
    {
        case VariablePositionType::Center:
            shift_x = 0.5;
            shift_y = 0.5;
            shift_z = 0.5;
            break;
        case VariablePositionType::XFace:
            shift_x = 0.0;
            shift_y = 0.5;
            shift_z = 0.5;
            break;
        case VariablePositionType::YFace:
            shift_x = 0.5;
            shift_y = 0.0;
            shift_z = 0.5;
            break;
        case VariablePositionType::ZFace:
            shift_x = 0.5;
            shift_y = 0.5;
            shift_z = 0.0;
            break;
        case VariablePositionType::Corner:
            shift_x = 0.0;
            shift_y = 0.0;
            shift_z = 0.0;
            break;
        default:
            // Or throw exception? Defaulting to Center for now or Null
            break;
    }

    for (auto& pair : field_map)
    {
        Domain3DUniform* s = pair.first;
        field3*          f = pair.second;

        OPENMP_PARALLEL_FOR()
        for (int i = 0; i < f->get_nx(); ++i)
        {
            for (int j = 0; j < f->get_ny(); ++j)
            {
                for (int k = 0; k < f->get_nz(); ++k)
                {
                    double gx     = s->get_offset_x() + (shift_x + i) * s->get_hx();
                    double gy     = s->get_offset_y() + (shift_y + j) * s->get_hy();
                    double gz     = s->get_offset_z() + (shift_z + k) * s->get_hz();
                    (*f)(i, j, k) = func(gx, gy, gz);
                }
            }
        }
    }
}

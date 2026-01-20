#include "variable.h"
#include <algorithm>
#include <array>

namespace
{
constexpr std::array<LocationType, 4> kBoundaryLocations2D = {
    LocationType::Left,
    LocationType::Right,
    LocationType::Down,
    LocationType::Up,
};
}

Variable::Variable(const std::string& in_name)
    : name(in_name)
{}

/**
 * @brief Attach this variable to a geometry.
 * @param g Geometry to bind to this variable.
 */
void Variable::set_geometry(Geometry2D& g)
{
    geometry = &g;
    for (auto& domainAdjPair : geometry->adjacency)
    {
        Domain2DUniform* domain = domainAdjPair.first;
        for (auto& locToNeighbor : domainAdjPair.second)
        {
            LocationType loc               = locToNeighbor.first;
            boundary_type_map[domain][loc] = PDEBoundaryType::Adjacented;
        }
    }

    // We suppose that when set the variable, the consstuction of the geometry has been finished
    // So we initial the corner boundary map here
    for (auto& domain : geometry->domains)
    {
        left_up_corner_value_map[domain]    = 0.0;
        right_down_corner_value_map[domain] = 0.0;
    }
}

/**
 * @brief Check geometry (before set fields)
 * @throws std::runtime_error If geometry is not set, domain is not part of the geometry,
 *         or the domain mesh size is invalid.
 */
void Variable::check_geometry(Domain2DUniform* s)
{
    if (geometry == nullptr)
        throw std::runtime_error("Variable has no geometry set");
    if (s->parent != geometry)
        throw std::runtime_error("Domain not found in geometry");
    if (s->nx <= 0 || s->ny <= 0)
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
void Variable::set_center_field(Domain2DUniform* s, field2& f)
{
    check_geometry(s);

    if (f.get_name() == "Default")
        f.init(s->nx, s->ny, name + "_" + s->name);
    else
        f.init(s->nx, s->ny);

    field_map[s] = &f;

    // Center variable only need left and down buffer
    buffer_map[s][LocationType::Left] = new double[s->ny];
    buffer_map[s][LocationType::Down] = new double[s->nx];

    position_type = VariablePositionType::Center;
}

void Variable::set_x_edge_field(Domain2DUniform* s, field2& f)
{
    check_geometry(s);

    if (f.get_name() == "Default")
        f.init(s->nx, s->ny, name + "_" + s->name);
    else
        f.init(s->nx, s->ny);

    field_map[s] = &f;

    buffer_map[s][LocationType::Left]  = new double[s->ny];
    buffer_map[s][LocationType::Right] = new double[s->ny];
    buffer_map[s][LocationType::Down]  = new double[s->nx];
    buffer_map[s][LocationType::Up]    = new double[s->nx];

    position_type = VariablePositionType::XFaceCenter;
}

void Variable::set_y_edge_field(Domain2DUniform* s, field2& f)
{
    check_geometry(s);

    if (f.get_name() == "Default")
        f.init(s->nx, s->ny, name + "_" + s->name);
    else
        f.init(s->nx, s->ny);

    field_map[s] = &f;

    buffer_map[s][LocationType::Left]  = new double[s->ny];
    buffer_map[s][LocationType::Right] = new double[s->ny];
    buffer_map[s][LocationType::Down]  = new double[s->nx];
    buffer_map[s][LocationType::Up]    = new double[s->nx];

    position_type = VariablePositionType::YFaceCenter;
}

void Variable::set_corner_field(Domain2DUniform* s, field2& f)
{
    check_geometry(s);

    if (f.get_name() == "Default")
        f.init(s->nx + 1, s->ny + 1, name + "_" + s->name);
    else
        f.init(s->nx + 1, s->ny + 1);

    field_map[s] = &f;

    buffer_map[s][LocationType::Left]  = new double[s->ny + 1];
    buffer_map[s][LocationType::Right] = new double[s->ny + 1];
    buffer_map[s][LocationType::Down]  = new double[s->nx + 1];
    buffer_map[s][LocationType::Up]    = new double[s->nx + 1];

    position_type = VariablePositionType::Corner;
}

void Variable::set_boundary_type(Domain2DUniform* s, LocationType loc, PDEBoundaryType type)
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

void Variable::set_boundary_type(Domain2DUniform*                                                s,
                                 std::initializer_list<std::pair<LocationType, PDEBoundaryType>> list)
{
    for (const auto& pair : list)
    {
        set_boundary_type(s, pair.first, pair.second);
    }
}

void Variable::fill_boundary_type(PDEBoundaryType type)
{
    if (geometry == nullptr)
        throw std::runtime_error("Variable has no geometry set");

    for (auto* domain : geometry->domains)
    {
        auto& type_map = boundary_type_map[domain];
        for (LocationType loc : kBoundaryLocations2D)
        {
            if (type_map.find(loc) != type_map.end())
                continue;
            set_boundary_type(domain, loc, type);
        }
    }
}

void Variable::set_boundary_value(Domain2DUniform* s, LocationType loc, double in_value)
{
    check_geometry(s);
    has_boundary_value_map[s][loc] = true;
    if (loc == LocationType::Left || loc == LocationType::Right)
    {
        boundary_value_map[s][loc] = new double[s->ny];
        for (int j = 0; j < s->ny; j++)
            boundary_value_map[s][loc][j] = in_value;

        if (loc == LocationType::Left)
            left_up_corner_value_map[s] = in_value;
    }
    else if (loc == LocationType::Down || loc == LocationType::Up)
    {
        boundary_value_map[s][loc] = new double[s->nx];
        for (int i = 0; i < s->nx; i++)
            boundary_value_map[s][loc][i] = in_value;

        if (loc == LocationType::Down)
            right_down_corner_value_map[s] = in_value;
    }
}

void Variable::set_boundary_value_from_func_global(Domain2DUniform*                      s,
                                                   LocationType                          loc,
                                                   std::function<double(double, double)> f)
{
    check_geometry(s);
    has_boundary_value_map[s][loc] = true;

    double shift_x = 0.5;
    double shift_y = 0.5;

    switch (position_type)
    {
        case VariablePositionType::Center:
            shift_x = 0.5;
            shift_y = 0.5;
            break;
        case VariablePositionType::XEdge:
            shift_x = 0.0;
            shift_y = 0.5;
            break;
        case VariablePositionType::YEdge:
            shift_x = 0.5;
            shift_y = 0.0;
            break;
        case VariablePositionType::Corner:
            shift_x = 0.0;
            shift_y = 0.0;
            break;
        default:
            break;
    }

    if (loc == LocationType::Left || loc == LocationType::Right)
    {
        int j_size;
        if (position_type == VariablePositionType::Corner)
            j_size = s->ny + 1;
        else
            j_size = s->ny;

        boundary_value_map[s][loc] = new double[j_size];

        // Determine i index (ghost node or boundary node)
        // Left: i = -1. Right: i = nx.
        int i_idx = (loc == LocationType::Left) ? -1 : s->nx;

        for (int j = 0; j < j_size; j++)
        {
            double gx                     = s->get_offset_x() + (shift_x + i_idx) * s->get_hx();
            double gy                     = s->get_offset_y() + (shift_y + j) * s->get_hy();
            boundary_value_map[s][loc][j] = f(gx, gy);
        }
    }
    else if (loc == LocationType::Down || loc == LocationType::Up)
    {
        int i_size;
        if (position_type == VariablePositionType::Corner)
            i_size = s->nx + 1;
        else
            i_size = s->nx;

        boundary_value_map[s][loc] = new double[i_size];

        // Down: j = -1. Up: j = ny.
        int j_idx = (loc == LocationType::Down) ? -1 : s->ny;

        for (int i = 0; i < i_size; i++)
        {
            double gx                     = s->get_offset_x() + (shift_x + i) * s->get_hx();
            double gy                     = s->get_offset_y() + (shift_y + j_idx) * s->get_hy();
            boundary_value_map[s][loc][i] = f(gx, gy);
        }
    }
}

void Variable::fill_boundary_value_from_func_global(std::function<double(double, double)> f)
{
    if (geometry == nullptr)
        throw std::runtime_error("Variable has no geometry set");

    for (auto* domain : geometry->domains)
    {
        auto& type_map = boundary_type_map[domain];
        auto& has_map  = has_boundary_value_map[domain];
        for (LocationType loc : kBoundaryLocations2D)
        {
            auto type_it = type_map.find(loc);
            if (type_it == type_map.end())
                continue;
            if (type_it->second == PDEBoundaryType::Adjacented || type_it->second == PDEBoundaryType::Null)
                continue;
            auto has_it = has_map.find(loc);
            if (has_it != has_map.end() && has_it->second)
                continue;

            set_boundary_value_from_func_global(domain, loc, f);
        }
    }
}

void Variable::set_value_from_func_global(std::function<double(double, double)> func)
{
    double shift_x = 0.5;
    double shift_y = 0.5;

    switch (position_type)
    {
        case VariablePositionType::Center:
            shift_x = 0.5;
            shift_y = 0.5;
            break;
        case VariablePositionType::XEdge:
            shift_x = 0.0;
            shift_y = 0.5;
            break;
        case VariablePositionType::YEdge:
            shift_x = 0.5;
            shift_y = 0.0;
            break;
        case VariablePositionType::Corner:
            shift_x = 0.0;
            shift_y = 0.0;
            break;
        default:
            // Or throw exception? Defaulting to Center for now or Null
            break;
    }

    for (auto& pair : field_map)
    {
        Domain2DUniform* s = pair.first;
        field2*          f = pair.second;

// Iterate over the field
#pragma omp parallel for collapse(2)
        for (int i = 0; i < f->get_nx(); ++i)
        {
            for (int j = 0; j < f->get_ny(); ++j)
            {
                double gx  = s->get_offset_x() + (shift_x + i) * s->get_hx();
                double gy  = s->get_offset_y() + (shift_y + j) * s->get_hy();
                (*f)(i, j) = func(gx, gy);
            }
        }
    }
}

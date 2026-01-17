#include "variable.h"
#include <algorithm>

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

    position_type = VariablePositionType::XEdge;
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

    position_type = VariablePositionType::YEdge;
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

void Variable::set_boundary_value(Domain2DUniform* s, LocationType loc, std::function<double(double)> f)
{
    check_geometry(s);
    // TODO:
}
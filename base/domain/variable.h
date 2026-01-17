#pragma once

#include "base/field/field2.h"
#include "domain2d.h"
#include "geometry2d.h"
#include <functional>
#include <stdexcept>
#include <unordered_map>

class Variable
{
public:
    std::string                                                                     name;
    Geometry2D*                                                                     geometry = nullptr;
    std::unordered_map<Domain2DUniform*, field2*>                                   field_map;
    std::unordered_map<Domain2DUniform*, std::unordered_map<LocationType, double*>> buffer_map;
    std::unordered_map<Domain2DUniform*, double> left_up_corner_value_map;    // Only for v, the value on the node
    std::unordered_map<Domain2DUniform*, double> right_down_corner_value_map; // Only for u，the value on the node

    std::unordered_map<Domain2DUniform*, std::unordered_map<LocationType, PDEBoundaryType>> boundary_type_map;
    std::unordered_map<Domain2DUniform*, std::unordered_map<LocationType, bool>>            has_boundary_value_map;
    std::unordered_map<Domain2DUniform*, std::unordered_map<LocationType, double*>>         boundary_value_map;

    Variable() = default;
    Variable(const std::string& in_name);
    ~Variable() = default;

    VariablePositionType position_type = VariablePositionType::Null;

    void set_geometry(Geometry2D& g);

    void check_geometry(Domain2DUniform* s);

    void set_center_field(Domain2DUniform* s, field2& f);
    void set_x_edge_field(Domain2DUniform* s, field2& f);
    void set_y_edge_field(Domain2DUniform* s, field2& f);
    void set_corner_field(Domain2DUniform* s, field2& f);

    void set_boundary_type(Domain2DUniform* s, LocationType loc, PDEBoundaryType type);
    void set_boundary_type(Domain2DUniform* s, std::initializer_list<std::pair<LocationType, PDEBoundaryType>> list);

    void set_boundary_value(Domain2DUniform* s, LocationType loc, double in_value);
    void
    set_boundary_value_from_func_global(Domain2DUniform* s, LocationType loc, std::function<double(double, double)> f);

    // void set_boundary_func_local(Domain2DUniform* s, LocationType loc, double* value);
    // void set_boundary_func_global(Domain2DUniform* s, LocationType loc, double* value);

    void set_value_from_func_global(std::function<double(double, double)> func);
};
#pragma once

#include <unordered_map>
#include <stdexcept>
#include "geometry2d.h"
#include "domain2d.h"
#include "base/field/field2.h"

class Variable
{
public:
    std::string  name;
    Geometry2D* geometry = nullptr;
    std::unordered_map<Domain2DUniform*, field2*> field_map;
    std::unordered_map<Domain2DUniform*, std::unordered_map<LocationType, field2*>> buffer_map;
    std::unordered_map<Domain2DUniform*, double> left_up_corner_map;    //Only for v
    std::unordered_map<Domain2DUniform*, double> right_down_corner_map; //Only for u

    std::unordered_map<Domain2DUniform*, std::unordered_map<LocationType, PDEBoundaryType>> boundary_type_map;

    Variable() = default;
    Variable(const std::string& in_name);
    ~Variable() = default;

    VariablePositionType position_type = VariablePositionType::Null;

    void set_geometry(Geometry2D& g);

    void check_geometry(Domain2DUniform* s);

    void set_center_field(Domain2DUniform* s, field2& f);
    void set_x_edge_field(Domain2DUniform* s, field2& f);
    void set_y_edge_field(Domain2DUniform* s, field2& f);

    void set_boundary_type(Domain2DUniform* s, LocationType loc, PDEBoundaryType type);
};
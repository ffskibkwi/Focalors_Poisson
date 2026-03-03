#pragma once

#include "base/field/field2.h"
#include "domain2d.h"
#include "geometry2d.h"
#include <functional>
#include <stdexcept>
#include <unordered_map>

class Variable2D
{
public:
    std::string                                                                     name;
    Geometry2D*                                                                     geometry = nullptr;
    std::unordered_map<Domain2DUniform*, field2*>                                   field_map;
    std::unordered_map<Domain2DUniform*, std::unordered_map<LocationType, double*>> buffer_map;
    std::unordered_map<Domain2DUniform*, double> xneg_ypos_corner_map; // Only for v, the value on the node
    std::unordered_map<Domain2DUniform*, double> xpos_yneg_corner_map; // Only for u，the value on the node

    std::unordered_map<Domain2DUniform*, std::unordered_map<LocationType, PDEBoundaryType>> boundary_type_map;
    std::unordered_map<Domain2DUniform*, std::unordered_map<LocationType, bool>>            has_boundary_value_map;
    std::unordered_map<Domain2DUniform*, std::unordered_map<LocationType, double*>>         boundary_value_map;

    Variable2D() = default;
    Variable2D(const std::string& in_name);
    virtual ~Variable2D();

    VariablePositionType position_type = VariablePositionType::Null;

    virtual void set_geometry(Geometry2D& g);

    virtual void check_geometry(Domain2DUniform* s);

    virtual void set_center_field(Domain2DUniform* s, field2& f);
    virtual void set_x_edge_field(Domain2DUniform* s, field2& f);
    virtual void set_y_edge_field(Domain2DUniform* s, field2& f);
    virtual void set_corner_field(Domain2DUniform* s, field2& f);
    virtual void set_inner_field(Domain2DUniform* s, field2& f);

    virtual void set_boundary_type(Domain2DUniform* s, LocationType loc, PDEBoundaryType type);
    virtual void set_boundary_type(Domain2DUniform*                                                s,
                                   std::initializer_list<std::pair<LocationType, PDEBoundaryType>> list);
    virtual void set_boundary_type(PDEBoundaryType type);

    virtual void set_boundary_value(Domain2DUniform* s, LocationType loc, double in_value);
    virtual void set_boundary_value(Domain2DUniform* s, LocationType loc, std::function<double(double, double)> f);
    virtual void set_boundary(std::function<double(double, double)> f);

    virtual void set_buffer_value(Domain2DUniform* s, LocationType loc, std::function<double(double, double)> f);
    virtual void set_corner(Domain2DUniform* s, std::function<double(double, double)> f);

    virtual void set_buffer(std::function<double(double, double)> f);
    virtual void set_corner(std::function<double(double, double)> f);

    virtual void set_value(std::function<double(double, double)> func);
};

#pragma once

#include "base/field/field3.h"
#include "domain3d.h"
#include "geometry3d.h"
#include <functional>
#include <stdexcept>
#include <unordered_map>

class Variable3D
{
public:
    std::string                                                                     name;
    Geometry3D*                                                                     geometry = nullptr;
    std::unordered_map<Domain3DUniform*, field3*>                                   field_map;
    std::unordered_map<Domain3DUniform*, std::unordered_map<LocationType, field2*>> buffer_map;
    std::unordered_map<Domain3DUniform*, double*> corner_x_map, corner_y_map, corner_z_map;

    std::unordered_map<Domain3DUniform*, std::unordered_map<LocationType, PDEBoundaryType>> boundary_type_map;
    std::unordered_map<Domain3DUniform*, std::unordered_map<LocationType, bool>>            has_boundary_value_map;
    std::unordered_map<Domain3DUniform*, std::unordered_map<LocationType, field2*>>         boundary_value_map;

    Variable3D() = default;
    Variable3D(const std::string& in_name);
    virtual ~Variable3D();

    VariablePositionType position_type = VariablePositionType::Null;

    virtual void set_geometry(Geometry3D& g);

    virtual void check_geometry(Domain3DUniform* s);

    virtual void set_center_field(Domain3DUniform* s, field3& f);
    virtual void set_x_face_center_field(Domain3DUniform* s, field3& f);
    virtual void set_y_face_center_field(Domain3DUniform* s, field3& f);
    virtual void set_z_face_center_field(Domain3DUniform* s, field3& f);

    virtual void set_boundary_type(Domain3DUniform* s, LocationType loc, PDEBoundaryType type);
    virtual void set_boundary_type(Domain3DUniform*                                                s,
                                   std::initializer_list<std::pair<LocationType, PDEBoundaryType>> list);
    virtual void set_boundary_type(PDEBoundaryType type);

    virtual void set_boundary_value(Domain3DUniform* s, LocationType loc, double in_value);
    virtual void
    set_boundary_value(Domain3DUniform* s, LocationType loc, std::function<double(double, double, double)> f);
    virtual void set_boundary(std::function<double(double, double, double)> f);

    virtual void
    set_buffer_value(Domain3DUniform* s, LocationType loc, std::function<double(double, double, double)> f);
    virtual void set_corner(Domain3DUniform* s, std::function<double(double, double, double)> f);

    virtual void set_buffer(std::function<double(double, double, double)> f);
    virtual void set_corner(std::function<double(double, double, double)> f);

    virtual void set_value(std::function<double(double, double, double)> func);
};

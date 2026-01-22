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
    std::unordered_map<Domain3DUniform*, double*> corner_value_map_x, corner_value_map_y, corner_value_map_z;

    std::unordered_map<Domain3DUniform*, std::unordered_map<LocationType, PDEBoundaryType>> boundary_type_map;
    std::unordered_map<Domain3DUniform*, std::unordered_map<LocationType, bool>>            has_boundary_value_map;
    std::unordered_map<Domain3DUniform*, std::unordered_map<LocationType, field2*>>         boundary_value_map;

    Variable3D() = default;
    Variable3D(const std::string& in_name);
    ~Variable3D();

    VariablePositionType position_type = VariablePositionType::Null;

    void set_geometry(Geometry3D& g);

    void check_geometry(Domain3DUniform* s);

    void set_center_field(Domain3DUniform* s, field3& f);
    void set_x_face_center_field(Domain3DUniform* s, field3& f);
    void set_y_face_center_field(Domain3DUniform* s, field3& f);
    void set_z_face_center_field(Domain3DUniform* s, field3& f);

    void set_boundary_type(Domain3DUniform* s, LocationType loc, PDEBoundaryType type);
    void set_boundary_type(Domain3DUniform* s, std::initializer_list<std::pair<LocationType, PDEBoundaryType>> list);
    void fill_boundary_type(PDEBoundaryType type);

    void set_boundary_value(Domain3DUniform* s, LocationType loc, double in_value);
    void set_boundary_value_from_func_global(Domain3DUniform*                              s,
                                             LocationType                                  loc,
                                             std::function<double(double, double, double)> f);
    void fill_boundary_value_from_func_global(std::function<double(double, double, double)> f);

    void set_value_from_func_global(std::function<double(double, double, double)> func);

private:
    void cleanup_buffers();
};

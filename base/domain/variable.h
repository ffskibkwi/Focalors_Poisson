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

    Variable() = default;
    Variable(const std::string& in_name);
    ~Variable() = default;

    void set_geometry(Geometry2D& g);

    void set_field(Domain2DUniform& s, field2& f);
};
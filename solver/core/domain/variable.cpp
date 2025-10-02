#include "variable.h"
#include <algorithm>

/**
 * @brief Attach this variable to a geometry.
 * @param g Geometry to bind to this variable.
 */
void Variable::set_geometry(Geometry2D& g)
{
    geometry = &g;
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
void Variable::set_field(Domain2DUniform& s, field2& f)
{
    if (geometry == nullptr)
        throw std::runtime_error("Variable has no geometry set");
    if (s.parent != geometry)
        throw std::runtime_error("Domain not found in geometry");
    if (s.nx <= 0 || s.ny <= 0)
        throw std::runtime_error("Domain mesh size invalid");

    if (f.get_name() == "Default")
        f.init(s.nx, s.ny, "variable");
    else
        f.init(s.nx, s.ny);

    field_map[&s] = &f;
}



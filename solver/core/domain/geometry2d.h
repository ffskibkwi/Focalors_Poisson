#pragma once

#include <vector>
#include <unordered_map>
#include <string>
#include <stdexcept>
#include "core/boundary/boundary_type.h"
#include "core/domain/domain2d.h"

class Geometry2D
{
public:
    // Subdomain storage and index
    std::vector<Domain2DUniform*> domains;
    std::unordered_map<Domain2DUniform*, std::unordered_map<LocationType, Domain2DUniform*>> adjacency;

    bool is_checked = false;
    Domain2DUniform* main_domain = nullptr;

    Geometry2D() = default;
    ~Geometry2D() = default;

    void add_domain(Domain2DUniform& s);
    void connect(Domain2DUniform& a, LocationType dir, Domain2DUniform& b);

    void check();

    void solve_prepare();
};
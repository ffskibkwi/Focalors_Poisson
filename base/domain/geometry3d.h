#pragma once

#include <vector>
#include <unordered_map>
#include <string>
#include <stdexcept>
#include <algorithm>
#include <queue>
#include <unordered_set>
#include "base/location_boundary.h"
#include "domain3d.h"
#include "geometry_tree.hpp"

class Geometry3D
{
public:
    // Subdomain storage and index
    std::vector<Domain3DUniform*> domains;
    std::unordered_map<Domain3DUniform*, std::unordered_map<LocationType, Domain3DUniform*>> adjacency;

    bool is_checked  = false;
    bool is_prepared = false;

    Domain3DUniform* tree_root = nullptr;
    std::unordered_map<Domain3DUniform*, std::unordered_map<LocationType, Domain3DUniform*>> tree_map;
    std::unordered_map<Domain3DUniform*, std::pair<LocationType, Domain3DUniform*>> parent_map;

    Geometry3D() = default;
    ~Geometry3D();

    void add_domain(Domain3DUniform* s);
    void connect(Domain3DUniform* a, LocationType dir, Domain3DUniform* b);

    void check();
    void solve_prepare();

    void set_position(Domain3DUniform* ref_domain, double pos_x, double pos_y, double pos_z);

private:
    // Check the single connectedness of the geometry
    bool is_single_connected() const;
    // Build the geometry tree
    void build_tree();
};

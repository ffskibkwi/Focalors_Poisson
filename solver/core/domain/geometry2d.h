#pragma once

#include <vector>
#include <unordered_map>
#include <string>
#include <stdexcept>
#include <algorithm>
#include <queue>
#include <unordered_set>
#include "core/boundary/boundary_type.h"
#include "domain2d.h"
#include "geometry_tree.hpp"

class Geometry2D
{
public:
    // Subdomain storage and index
    std::vector<Domain2DUniform*> domains;
    std::unordered_map<Domain2DUniform*, std::unordered_map<LocationType, Domain2DUniform*>> adjacency;

    bool is_checked = false;
    // Domain2DUniform* main_domain = nullptr; // deprecated: tree-based analysis replaces single main domain

    // using GeometryTreeNode2D = GeometryTreeNode<Domain2DUniform> (in geometry_tree.hpp)
    // GeometryTreeNode2D* tree_root = nullptr; (for geometry_tree_old.hpp)

    Domain2DUniform* tree_root = nullptr;
    std::unordered_map<Domain2DUniform*, std::vector<Domain2DUniform*>> tree_map;

    Geometry2D() = default;
    ~Geometry2D();

    void add_domain(Domain2DUniform& s);
    void connect(Domain2DUniform& a, LocationType dir, Domain2DUniform& b);

    void check();

    void solve_prepare();

private:
    // Check the single connectedness of the geometry
    bool is_single_connected() const;
    // Build the geometry tree
    void build_tree();
};
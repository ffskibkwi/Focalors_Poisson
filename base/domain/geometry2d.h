#pragma once

#include "base/location_boundary.h"
#include "domain2d.h"
#include "geometry_tree.hpp"
#include <algorithm>
#include <initializer_list>
#include <queue>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

class Geometry2D
{
public:
    // Subdomain storage and index
    std::vector<Domain2DUniform*>                                                            domains;
    std::unordered_map<Domain2DUniform*, std::unordered_map<LocationType, Domain2DUniform*>> adjacency;

    bool is_checked  = false;
    bool is_prepared = false;
    // Domain2DUniform* main_domain = nullptr; // deprecated: tree-based analysis replaces single main domain

    // using GeometryTreeNode2D = GeometryTreeNode<Domain2DUniform> (in geometry_tree.hpp)
    // GeometryTreeNode2D* tree_root = nullptr; (for geometry_tree_old.hpp)

    Domain2DUniform*                                                                         tree_root = nullptr;
    std::unordered_map<Domain2DUniform*, std::unordered_map<LocationType, Domain2DUniform*>> tree_map;
    std::unordered_map<Domain2DUniform*, std::pair<LocationType, Domain2DUniform*>>          parent_map;

    Geometry2D() = default;
    ~Geometry2D();

    void add_domain(Domain2DUniform* s);
    void add_domain(std::initializer_list<Domain2DUniform*> list);
    void connect(Domain2DUniform& a, LocationType dir, Domain2DUniform& b);

    void check();
    void solve_prepare();

    void set_global_spatial_step(double hx, double hy);

    void axis(Domain2DUniform* d, LocationType loc);

    void update_offset_x(Domain2DUniform* d);
    void update_offset_y(Domain2DUniform* d);

    void global_move_x(double x);
    void global_move_y(double y);

    void set_position(Domain2DUniform* ref_domain, double pos_x, double pos_y);

private:
    // Check the single connectedness of the geometry
    bool is_single_connected() const;
    // Build the geometry tree
    void build_tree();
};
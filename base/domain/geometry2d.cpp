#include "geometry2d.h"

Geometry2D::~Geometry2D() {}

/**
 * @brief Add a domain to the geometry if not already present.
 * @param s Domain to add.
 *
 * Ensures uniqueness of `s` within the geometry and sets its `parent`
 * pointer to this `Geometry2D` instance.
 */
void Geometry2D::add_domain(Domain2DUniform* s)
{
    if (s == nullptr)
    {
        throw std::invalid_argument("add_domain: nullptr");
        return;
    }
    if (std::find(domains.begin(), domains.end(), s) == domains.end())
        domains.push_back(s);
    s->parent = this;
}

void Geometry2D::add_domain(std::initializer_list<Domain2DUniform*> list)
{
    for (auto* s : list)
    {
        if (s == nullptr)
        {
            throw std::invalid_argument("add_domain: nullptr");
            return;
        }
        add_domain(s);
    }
}

/**
 * @brief Connect two domains with a directional adjacency.
 * @param a Source domain.
 * @param dir Direction from domain `a` to domain `b`.
 * @param b Target domain.
 * @throws std::invalid_argument If `dir` is `Front` or `Back` which are unsupported in 2D.
 *
 * The connection is stored symmetrically: `a --dir--> b` and `b --opposite(dir)--> a`.
 * Default boundary types on the connected faces are set to Dirichlet for both domains.
 */
void Geometry2D::connect(Domain2DUniform* a, LocationType dir, Domain2DUniform* b)
{
    // Here defaultly add the first domain and the second domain into geo
    // Maybe it is dangerous~
    add_domain({a, b});
    if (dir == LocationType::Front || dir == LocationType::Back)
        throw std::invalid_argument("Geometry2D does not support Front/Back");

    adjacency[a][dir]           = b;
    adjacency[b][opposite(dir)] = a;

    // 不再在几何层面设置边界类型，边界类型应绑定在 Variable 上

    // Set the size
    // In  this function, b is decided by a
    if (dir == LocationType::Left || dir == LocationType::Right)
    {
        b->set_ly(a->ly);
        b->set_ny(a->ny);
    }
    else if (dir == LocationType::Up || dir == LocationType::Down)
    {
        b->set_lx(a->lx);
        b->set_nx(a->nx);
    }
}

/**
 * @brief Validate profiles and boundaries, and determine the main domain.
 * @throws std::runtime_error If a domain has invalid profile/boundary,
 *         if multiple main domains are detected, or if none is found.
 */
void Geometry2D::check()
{
    for (auto* s : domains)
    {
        // Check the domain profile
        if (!s->check_profile())
            throw std::runtime_error("Domain " + s->name + " has invalid profile");
        // 不再在几何层面校验边界类型（边界由 Variable 管理）
        // Deprecated: single main domain detection is removed in favor of tree-based analysis
    }

    // Check the single connectedness of the geometry
    if (!is_single_connected())
        throw std::runtime_error("Geometry is not single connected");
    is_checked = true;
}

/**
 * @brief Prepare geometry for solving after validation.
 * @throws std::runtime_error If `check()` has not been called.
 */
void Geometry2D::solve_prepare()
{
    if (!is_checked)
        throw std::runtime_error("Geometry2D is not checked");

    // Build the geometry tree, for later use
    if (!is_prepared)
        build_tree();
}
bool Geometry2D::is_single_connected() const
{
    if (domains.empty())
        return true;
    // Start from any node in adjacency; if adjacency is empty but domains>1, it is not connected
    Domain2DUniform* start = nullptr;
    if (!adjacency.empty())
        start = adjacency.begin()->first;
    else
        return domains.size() == 1;

    std::unordered_set<Domain2DUniform*> visited;
    std::queue<Domain2DUniform*>         q;
    q.push(start);
    visited.insert(start);

    auto enqueue = [&](Domain2DUniform* a, Domain2DUniform* b) {
        if (b && !visited.count(b))
        {
            visited.insert(b);
            q.push(b);
        }
    };

    while (!q.empty())
    {
        Domain2DUniform* u = q.front();
        q.pop();
        // Forward adjacency
        auto it = adjacency.find(u);
        if (it != adjacency.end())
        {
            for (const auto& kv : it->second)
                enqueue(u, kv.second);
        }
        // Reverse adjacency: traverse other keys, if its adjacency value contains u, it is also considered as an
        // undirected connection
        for (const auto& ap : adjacency)
        {
            for (const auto& kv : ap.second)
            {
                if (kv.second == u)
                    enqueue(ap.first, ap.first);
            }
        }
    }

    // All domains must be visited
    for (auto* s : domains)
    {
        if (!visited.count(s))
            return false;
    }
    return true;
}

void Geometry2D::build_tree()
{
    // for geometry_tree_old.hpp
    // TreeBuilder2D builder;
    // tree_root = builder.buildOptimalTree(adjacency);

    tree_root  = TreeUtils::findOptimalRoot(adjacency);
    tree_map   = TreeUtils::buildTreeMapFromRoot(tree_root, adjacency);
    parent_map = TreeUtils::buildParentMapFromTree(tree_map);

    is_prepared = true;
}

void Geometry2D::set_global_spatial_step(double hx, double hy)
{
    for (auto* s : domains)
    {
        s->set_spatial_step(hx, hy);
    }
}

void Geometry2D::set_position(Domain2DUniform* ref_domain, double pos_x, double pos_y)
{
    std::queue<Domain2DUniform*> q;

    q.push(ref_domain);
    while (!q.empty())
    {
        Domain2DUniform* currentDomain = q.front();
        q.pop();
        currentDomain->set_position(pos_x, pos_y);
        //
    }
}

void Geometry2D::axis(Domain2DUniform* d, LocationType loc)
{
    if (std::find(domains.begin(), domains.end(), d) == domains.end())
        throw std::runtime_error("axis: Domain is not in geometry");

    // Check direction type to decide if we are setting X or Y
    bool is_x_axis = (loc == LocationType::Left || loc == LocationType::Right);

    // Initialize offset for the starting domain
    if (is_x_axis)
    {
        if (loc == LocationType::Left)
            d->set_offset_x(0.0);
        else // Right
            d->set_offset_x(-d->get_lx());

        update_offset_x(d);
    }
    else
    {
        if (loc == LocationType::Down)
            d->set_offset_y(0.0);
        else if (loc == LocationType::Up)
            d->set_offset_y(-d->get_ly());
        else
            throw std::runtime_error("axis: Invalid LocationType for 2D");

        update_offset_y(d);
    }
}

void Geometry2D::update_offset_x(Domain2DUniform* d)
{
    // BFS to propagate offsets
    std::unordered_set<Domain2DUniform*> visited;
    std::queue<Domain2DUniform*>         q;
    q.push(d);
    visited.insert(d);

    while (!q.empty())
    {
        Domain2DUniform* u = q.front();
        q.pop();

        auto it = adjacency.find(u);
        if (it == adjacency.end())
            continue;

        for (const auto& pair : it->second)
        {
            LocationType     dir = pair.first;
            Domain2DUniform* v   = pair.second;

            if (visited.count(v))
                continue;

            double current_x = u->get_offset_x();
            double new_x     = current_x;
            if (dir == LocationType::Right)
                new_x = current_x + u->get_lx();
            else if (dir == LocationType::Left)
                new_x = current_x - v->get_lx();
            // For Up/Down, X offset is propagated unchanged (assuming alignment)
            v->set_offset_x(new_x);

            visited.insert(v);
            q.push(v);
        }
    }
}

void Geometry2D::update_offset_y(Domain2DUniform* d)
{
    // BFS to propagate offsets
    std::unordered_set<Domain2DUniform*> visited;
    std::queue<Domain2DUniform*>         q;
    q.push(d);
    visited.insert(d);

    while (!q.empty())
    {
        Domain2DUniform* u = q.front();
        q.pop();

        auto it = adjacency.find(u);
        if (it == adjacency.end())
            continue;

        for (const auto& pair : it->second)
        {
            LocationType     dir = pair.first;
            Domain2DUniform* v   = pair.second;

            if (visited.count(v))
                continue;

            double current_y = u->get_offset_y();
            double new_y     = current_y;
            if (dir == LocationType::Up)
                new_y = current_y + u->get_ly();
            else if (dir == LocationType::Down)
                new_y = current_y - v->get_ly();
            // For Left/Right, Y offset is propagated unchanged
            v->set_offset_y(new_y);

            visited.insert(v);
            q.push(v);
        }
    }
}

void Geometry2D::global_move_x(double x)
{
    for (auto* d : domains)
    {
        d->set_offset_x(d->get_offset_x() + x);
    }
}

void Geometry2D::global_move_y(double y)
{
    for (auto* d : domains)
    {
        d->set_offset_y(d->get_offset_y() + y);
    }
}
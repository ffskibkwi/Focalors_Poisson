#include "geometry3d.h"

Geometry3D::~Geometry3D() {}

/**
 * @brief Add a domain to the geometry if not already present.
 * @param s Domain to add.
 *
 * Ensures uniqueness of `s` within the geometry and sets its `parent`
 * pointer to this `Geometry3D` instance.
 */
void Geometry3D::add_domain(Domain3DUniform* s)
{
    if (std::find(domains.begin(), domains.end(), s) == domains.end())
        domains.push_back(s);
    s->parent = this;
}

/**
 * @brief Connect two domains with a directional adjacency.
 * @param a Source domain.
 * @param dir Direction from domain `a` to domain `b`.
 * @param b Target domain.
 *
 * The connection is stored symmetrically: `a --dir--> b` and `b --opposite(dir)--> a`.
 */
void Geometry3D::connect(Domain3DUniform* a, LocationType dir, Domain3DUniform* b)
{
    // Here defaultly add the first domain and the second domain into geo
    add_domain(a);
    add_domain(b);

    adjacency[a][dir]           = b;
    adjacency[b][opposite(dir)] = a;

    // Set the size
    // In this function, b is decided by a
    if (dir == LocationType::Left || dir == LocationType::Right)
    {
        b->set_ly(a->ly);
        b->set_ny(a->ny);
        b->set_lz(a->lz);
        b->set_nz(a->nz);
    }
    else if (dir == LocationType::Front || dir == LocationType::Back)
    {
        b->set_lx(a->lx);
        b->set_nx(a->nx);
        b->set_lz(a->lz);
        b->set_nz(a->nz);
    }
    else if (dir == LocationType::Down || dir == LocationType::Up)
    {
        b->set_lx(a->lx);
        b->set_nx(a->nx);
        b->set_ly(a->ly);
        b->set_ny(a->ny);
    }
}

/**
 * @brief Validate profiles and boundaries, and determine the main domain.
 * @throws std::runtime_error If a domain has invalid profile/boundary,
 *         if multiple main domains are detected, or if none is found.
 */
void Geometry3D::check()
{
    for (auto* s : domains)
    {
        // Check the domain profile
        if (!s->check_profile())
            throw std::runtime_error("Domain " + s->name + " has invalid profile");
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
void Geometry3D::solve_prepare()
{
    if (!is_checked)
        throw std::runtime_error("Geometry3D is not checked");

    // Build the geometry tree, for later use
    if (!is_prepared)
        build_tree();
}

bool Geometry3D::is_single_connected() const
{
    if (domains.empty())
        return true;
    // Start from any node in adjacency; if adjacency is empty but domains>1, it is not connected
    Domain3DUniform* start = nullptr;
    if (!adjacency.empty())
        start = adjacency.begin()->first;
    else
        return domains.size() == 1;

    std::unordered_set<Domain3DUniform*> visited;
    std::queue<Domain3DUniform*>         q;
    q.push(start);
    visited.insert(start);

    auto enqueue = [&](Domain3DUniform* a, Domain3DUniform* b) {
        if (b && !visited.count(b))
        {
            visited.insert(b);
            q.push(b);
        }
    };

    while (!q.empty())
    {
        Domain3DUniform* u = q.front();
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

void Geometry3D::build_tree()
{
    if (domains.size() == 1)
        tree_root = domains[0];
    else
        tree_root = TreeUtils::findOptimalRoot(adjacency);
    tree_map                  = TreeUtils::buildTreeMapFromRoot(tree_root, adjacency);
    parent_map                = TreeUtils::buildParentMapFromTree(tree_map);
    hierarchical_solve_levels = TreeUtils::BFS(tree_root, tree_map);

    is_prepared = true;
}

void Geometry3D::set_position(Domain3DUniform* ref_domain, double pos_x, double pos_y, double pos_z)
{
    std::queue<Domain3DUniform*> q;

    q.push(ref_domain);
    while (!q.empty())
    {
        Domain3DUniform* currentDomain = q.front();
        q.pop();
        currentDomain->set_position(pos_x, pos_y, pos_z);
        // TODO: Update position based on adjacency
    }
}

void Geometry3D::set_global_spatial_step(double hx, double hy, double hz)
{
    for (auto* s : domains)
    {
        s->set_spatial_step(hx, hy, hz);
    }
}

void Geometry3D::axis(Domain3DUniform* d, LocationType loc)
{
    if (std::find(domains.begin(), domains.end(), d) == domains.end())
        throw std::runtime_error("axis: Domain is not in geometry");

    // Check direction type to decide if we are setting X or Y
    bool is_x_axis = (loc == LocationType::Left || loc == LocationType::Right);
    bool is_y_axis = (loc == LocationType::Front || loc == LocationType::Back);

    // Initialize offset for the starting domain
    if (is_x_axis)
    {
        if (loc == LocationType::Left)
            d->set_offset_x(0.0);
        else // Right
            d->set_offset_x(-d->get_lx());

        update_offset_x(d);
    }
    else if (is_y_axis)
    {
        if (loc == LocationType::Front)
            d->set_offset_y(0.0);
        else if (loc == LocationType::Back)
            d->set_offset_y(-d->get_ly());

        update_offset_y(d);
    }
    else
    {
        if (loc == LocationType::Down)
            d->set_offset_z(0.0);
        else if (loc == LocationType::Up)
            d->set_offset_z(-d->get_lz());

        update_offset_z(d);
    }
}

void Geometry3D::update_offset_x(Domain3DUniform* d)
{
    // BFS to propagate offsets
    std::unordered_set<Domain3DUniform*> visited;
    std::queue<Domain3DUniform*>         q;
    q.push(d);
    visited.insert(d);

    while (!q.empty())
    {
        Domain3DUniform* u = q.front();
        q.pop();

        auto it = adjacency.find(u);
        if (it == adjacency.end())
            continue;

        for (const auto& pair : it->second)
        {
            LocationType     dir = pair.first;
            Domain3DUniform* v   = pair.second;

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

void Geometry3D::update_offset_y(Domain3DUniform* d)
{
    // BFS to propagate offsets
    std::unordered_set<Domain3DUniform*> visited;
    std::queue<Domain3DUniform*>         q;
    q.push(d);
    visited.insert(d);

    while (!q.empty())
    {
        Domain3DUniform* u = q.front();
        q.pop();

        auto it = adjacency.find(u);
        if (it == adjacency.end())
            continue;

        for (const auto& pair : it->second)
        {
            LocationType     dir = pair.first;
            Domain3DUniform* v   = pair.second;

            if (visited.count(v))
                continue;

            double current_y = u->get_offset_y();
            double new_y     = current_y;
            if (dir == LocationType::Back)
                new_y = current_y + u->get_ly();
            else if (dir == LocationType::Front)
                new_y = current_y - v->get_ly();
            // For Left/Right, Y offset is propagated unchanged
            v->set_offset_y(new_y);

            visited.insert(v);
            q.push(v);
        }
    }
}

void Geometry3D::update_offset_z(Domain3DUniform* d)
{
    // BFS to propagate offsets
    std::unordered_set<Domain3DUniform*> visited;
    std::queue<Domain3DUniform*>         q;
    q.push(d);
    visited.insert(d);

    while (!q.empty())
    {
        Domain3DUniform* u = q.front();
        q.pop();

        auto it = adjacency.find(u);
        if (it == adjacency.end())
            continue;

        for (const auto& pair : it->second)
        {
            LocationType     dir = pair.first;
            Domain3DUniform* v   = pair.second;

            if (visited.count(v))
                continue;

            double current_z = u->get_offset_z();
            double new_z     = current_z;
            if (dir == LocationType::Up)
                new_z = current_z + u->get_lz();
            else if (dir == LocationType::Down)
                new_z = current_z - v->get_lz();
            // For Left/Right, Y offset is propagated unchanged
            v->set_offset_z(new_z);

            visited.insert(v);
            q.push(v);
        }
    }
}

void Geometry3D::global_move_x(double x)
{
    for (auto* d : domains)
    {
        d->set_offset_x(d->get_offset_x() + x);
    }
}

void Geometry3D::global_move_y(double y)
{
    for (auto* d : domains)
    {
        d->set_offset_y(d->get_offset_y() + y);
    }
}

void Geometry3D::global_move_z(double z)
{
    for (auto* d : domains)
    {
        d->set_offset_z(d->get_offset_z() + z);
    }
}

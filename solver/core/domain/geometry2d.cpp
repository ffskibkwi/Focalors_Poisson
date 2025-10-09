#include "geometry2d.h"
#include <algorithm>
#include <queue>
#include <unordered_set>
Geometry2D::~Geometry2D()
{
    if (tree_root) { deleteTree(tree_root); tree_root = nullptr; }
}


/**
 * @brief Add a domain to the geometry if not already present.
 * @param s Domain to add.
 *
 * Ensures uniqueness of `s` within the geometry and sets its `parent`
 * pointer to this `Geometry2D` instance.
 */
void Geometry2D::add_domain(Domain2DUniform& s)
{
    if (std::find(domains.begin(), domains.end(), &s) == domains.end())
        domains.push_back(&s);
    s.parent = this;
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
void Geometry2D::connect(Domain2DUniform& a, LocationType dir, Domain2DUniform& b)
{
    //Here defaultly add the first domain and the second domain into geo
    //Maybe it is dangerous~
    add_domain(a);
    add_domain(b);
    if (dir == LocationType::Front || dir == LocationType::Back)
        throw std::invalid_argument("Geometry2D does not support Front/Back");

    adjacency[&a][dir] = &b;
    adjacency[&b][opposite(dir)] = &a;

    //Set the boundary of the two domains
    a.set_boundary(dir, PDEBoundaryType::Dirichlet);
    b.set_boundary(opposite(dir), PDEBoundaryType::Dirichlet);

    //Set the size
    //In  this function, b is decided by a
    if (dir == LocationType::Left || dir == LocationType::Right)
    {
        b.set_ly(a.ly);
        b.set_ny(a.ny);
    }else if (dir == LocationType::Up || dir == LocationType::Down)
    {
        b.set_lx(a.lx);
        b.set_nx(a.nx);
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
        //Check the domain profile and boundary
        if (!s->check_profile())
            throw std::runtime_error("Domain " + s->name + " has invalid profile");
        if (!s->check_boundary())
            throw std::runtime_error("Domain " + s->name + " has invalid boundary");
        // Deprecated: single main domain detection is removed in favor of tree-based analysis
    }

    //Check the single connectedness of the geometry
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

    //Build the geometry tree, for later use
    build_tree();
}
bool Geometry2D::is_single_connected() const
{
    if (domains.empty()) return true;
    //Start from any node in adjacency; if adjacency is empty but domains>1, it is not connected
    Domain2DUniform* start = nullptr;
    if (!adjacency.empty()) start = adjacency.begin()->first;
    else return domains.size() == 1;

    std::unordered_set<Domain2DUniform*> visited;
    std::queue<Domain2DUniform*> q;
    q.push(start);
    visited.insert(start);

    auto enqueue = [&](Domain2DUniform* a, Domain2DUniform* b){
        if (b && !visited.count(b)) { visited.insert(b); q.push(b); }
    };

    while (!q.empty()) {
        Domain2DUniform* u = q.front(); q.pop();
        //Forward adjacency
        auto it = adjacency.find(u);
        if (it != adjacency.end()) {
            for (const auto& kv : it->second) enqueue(u, kv.second);
        }
        //Reverse adjacency: traverse other keys, if its adjacency value contains u, it is also considered as an undirected connection
        for (const auto& ap : adjacency) {
            for (const auto& kv : ap.second) {
                if (kv.second == u) enqueue(ap.first, ap.first);
            }
        }
    }

    //All domains must be visited
    for (auto* s : domains) {
        if (!visited.count(s)) return false;
    }
    return true;
}

void Geometry2D::build_tree()
{
    TreeBuilder2D builder;
    tree_root = builder.buildOptimalTree(adjacency);
}



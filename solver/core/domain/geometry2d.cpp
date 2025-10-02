#include "geometry2d.h"
#include <algorithm>

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
        
        //Check the domain degree and find the main domain
        size_t deg = 0;
        auto it = adjacency.find(s);
        if (it != adjacency.end()) deg = it->second.size();
        if (deg > 1)
        {
            if (main_domain != nullptr)
                throw std::runtime_error("Multiple main domains are not supported");
            main_domain = s;
        }
    }
    if (main_domain == nullptr)
        throw std::runtime_error("No main domain found");
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
}



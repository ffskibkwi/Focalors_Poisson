#pragma once

#include "d_domain_2d.h"
#include "base/field/field2.h"
#include <string>

/**
 * @brief Distributed Field 2D
 * 
 * Represents a scalar field distributed across the DDomain2D.
 * Holds the local slab of data and provides methods for scatter/gather.
 */
class DField2D
{
public:
    DField2D(DDomain2D* in_domain, std::string name = "d_field");
    ~DField2D() = default;

    // Access to local data
    field2& get_local_data() { return local_data; }
    const field2& get_local_data() const { return local_data; }
    
    DDomain2D* get_domain() const { return domain; }

    // Initialize local data to zero
    void zeros();

    // Scatter from a global field on root (rank 0) to distributed local fields
    // global_field is only significant on root
    void scatter_from_global(const field2& global_field);

    // Gather distributed local fields to a global field on root
    // global_field is resized on root to hold the result
    void gather_to_global(field2& global_field) const;

    // Helper to get global coordinate of a local point (i_local, j_local)
    // i_local: [0, local_nx), j_local: [0, local_ny)
    // Returns global (i, j)
    void get_global_coord(int i_local, int j_local, int& i_global, int& j_global) const;

private:
    DDomain2D* domain = nullptr;
    field2 local_data;
    std::string name;
};

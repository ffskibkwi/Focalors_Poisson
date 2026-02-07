#pragma once

#include "domain2d_mpi.h"
#include "variable2d.h"

#include <mpi.h>

class Variable2DSlabX : public Variable2D
{
public:
    Variable2DSlabX() = default;
    Variable2DSlabX(const std::string& in_name, MPI_Comm _communicator = MPI_COMM_WORLD);
    ~Variable2DSlabX();

    void set_geometry(Geometry2D& g);

    void set_center_field(Domain2DUniformMPI* s, field2& f);
    void set_x_edge_field(Domain2DUniformMPI* s, field2& f);
    void set_y_edge_field(Domain2DUniformMPI* s, field2& f);
    void set_corner_field(Domain2DUniformMPI* s, field2& f);

    void set_boundary_value(Domain2DUniformMPI* s, LocationType loc, double in_value);
    void set_boundary_value_from_func_global(Domain2DUniformMPI*                   s,
                                             LocationType                          loc,
                                             std::function<double(double, double)> f);
    void fill_boundary_value_from_func_global(std::function<double(double, double)> f);

    void set_value_from_func_global(std::function<double(double, double)> func);

    void print_slab_info();

    // slab
    MPI_Comm communicator = MPI_COMM_WORLD;
    int      mpi_rank, mpi_size;

    // hierarchical var contains tree root

    std::vector<Domain2DUniform*> hierarchical_slab_parents;
    std::vector<MPI_Comm>         hierarchical_slab_comms;
    std::vector<int>              hierarchical_slab_ranks;
    std::vector<int>              hierarchical_slab_sizes;
    std::vector<int>              hierarchical_slab_nxs;
    std::vector<int>              hierarchical_slab_disps;

    std::unordered_map<int, int> slab_parent_to_level;
};

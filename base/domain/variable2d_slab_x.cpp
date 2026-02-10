#include "variable2d_slab_x.h"
#include "base/parallel/mpi/mpi_misc.h"
#include "base/parallel/omp/enable_openmp.h"

#include <algorithm>
#include <array>

Variable2DSlabX::Variable2DSlabX(const std::string& in_name, MPI_Comm _communicator)
    : Variable2D(in_name)
    , communicator(_communicator)
{}

Variable2DSlabX::~Variable2DSlabX()
{
    // Usually we finalize mpi before we exit the process
    // for (auto comm : hierarchical_slab_comms)
    //     MPI_Comm_free(&comm);
}

/**
 * @brief Attach this variable to a geometry.
 * @param g Geometry to bind to this variable.
 */
void Variable2DSlabX::set_geometry(Geometry2D& g)
{
    MPI_Comm_rank(communicator, &mpi_rank);
    MPI_Comm_size(communicator, &mpi_size);

    Variable2D::set_geometry(g);

    for (auto level : g.hierarchical_solve_levels)
    {
        int nx_sum = 0;
        for (auto domain : level)
            nx_sum += domain->get_nx();

        int color = 0;
        int right = 0;
        for (int i = 0; i < level.size() - 1; i++)
        {
            auto domain = level[i];

            int comm_size = (double)domain->get_nx() / nx_sum * mpi_size;
            if (comm_size == 0)
                comm_size = 1;
            right += comm_size;

            if (mpi_rank >= right)
                color++;
        }

        MPI_Comm slab_comm;
        MPI_Comm_split(MPI_COMM_WORLD, color, mpi_rank, &slab_comm);

        auto current_domain = level[color];
        hierarchical_slab_parents.push_back(current_domain);
        hierarchical_slab_comms.push_back(slab_comm);

        int hierarchical_rank = 0;
        int hierarchical_size = 0;
        MPI_Comm_rank(slab_comm, &hierarchical_rank);
        MPI_Comm_size(slab_comm, &hierarchical_size);
        hierarchical_slab_ranks.push_back(hierarchical_rank);
        hierarchical_slab_sizes.push_back(hierarchical_size);
        hierarchical_slab_nxs.push_back(
            MPIUtils::get_slab_length(current_domain->nx, hierarchical_rank, hierarchical_size));
        hierarchical_slab_disps.push_back(
            MPIUtils::get_slab_displacement(current_domain->nx, hierarchical_rank, hierarchical_size));
    }

    // for O(1) access
    for (int i = 0; i < hierarchical_slab_parents.size(); i++)
    {
        auto domain_mpi                              = static_cast<Domain2DUniformMPI*>(hierarchical_slab_parents[i]);
        slab_parent_to_level[domain_mpi->get_uuid()] = i;
    }
}

/**
 * @brief Initialize a field and bind it to a domain within the attached geometry.
 * @param s Domain to which the field belongs.
 * @param f Field storage to initialize and register.
 * @throws std::runtime_error If geometry is not set, domain is not part of the geometry,
 *         or the domain mesh size is invalid.
 *
 * If the field name is "Default", it will be initialized with the name "variable".
 */
void Variable2DSlabX::set_center_field(Domain2DUniformMPI* s, field2& f)
{
    check_geometry(s);

    int level = 0;
    if (slab_parent_to_level.find(s->get_uuid()) != slab_parent_to_level.end())
        level = slab_parent_to_level[s->get_uuid()];
    else
        return;

    int snx = hierarchical_slab_nxs[level];
    int ny  = s->ny;

    f.init(snx, ny, name + "_" + s->name);

    field_map[s] = &f;

    // Center variable only need left and down buffer
    buffer_map[s][LocationType::Left] = new double[ny];
    buffer_map[s][LocationType::Down] = new double[snx];

    position_type = VariablePositionType::Center;
}

void Variable2DSlabX::set_x_edge_field(Domain2DUniformMPI* s, field2& f)
{
    check_geometry(s);

    int level = 0;
    if (slab_parent_to_level.find(s->get_uuid()) != slab_parent_to_level.end())
        level = slab_parent_to_level[s->get_uuid()];
    else
        return;

    int snx = hierarchical_slab_nxs[level];
    int ny  = s->ny;

    f.init(snx, ny, name + "_" + s->name);

    field_map[s] = &f;

    buffer_map[s][LocationType::Left]  = new double[ny];
    buffer_map[s][LocationType::Right] = new double[ny];
    buffer_map[s][LocationType::Down]  = new double[snx];
    buffer_map[s][LocationType::Up]    = new double[snx];

    position_type = VariablePositionType::XFace;
}

void Variable2DSlabX::set_y_edge_field(Domain2DUniformMPI* s, field2& f)
{
    check_geometry(s);

    int level = 0;
    if (slab_parent_to_level.find(s->get_uuid()) != slab_parent_to_level.end())
        level = slab_parent_to_level[s->get_uuid()];
    else
        return;

    int snx = hierarchical_slab_nxs[level];
    int ny  = s->ny;

    f.init(snx, ny, name + "_" + s->name);

    field_map[s] = &f;

    buffer_map[s][LocationType::Left]  = new double[ny];
    buffer_map[s][LocationType::Right] = new double[ny];
    buffer_map[s][LocationType::Down]  = new double[snx];
    buffer_map[s][LocationType::Up]    = new double[snx];

    position_type = VariablePositionType::YFace;
}

void Variable2DSlabX::set_corner_field(Domain2DUniformMPI* s, field2& f)
{
    check_geometry(s);

    int level = 0;
    if (slab_parent_to_level.find(s->get_uuid()) != slab_parent_to_level.end())
        level = slab_parent_to_level[s->get_uuid()];
    else
        return;

    int snx = hierarchical_slab_nxs[level];
    int ny  = s->ny + 1;

    f.init(snx, ny, name + "_" + s->name);

    field_map[s] = &f;

    buffer_map[s][LocationType::Left]  = new double[ny];
    buffer_map[s][LocationType::Right] = new double[ny];
    buffer_map[s][LocationType::Down]  = new double[snx];
    buffer_map[s][LocationType::Up]    = new double[snx];

    position_type = VariablePositionType::Corner;
}

void Variable2DSlabX::set_inner_field(Domain2DUniformMPI* s, field2& f)
{
    check_geometry(s);

    int level = 0;
    if (slab_parent_to_level.find(s->get_uuid()) != slab_parent_to_level.end())
        level = slab_parent_to_level[s->get_uuid()];
    else
        return;

    int snx = hierarchical_slab_nxs[level];
    int ny  = s->ny;

    f.init(snx, ny, name + "_" + s->name);

    field_map[s] = &f;

    buffer_map[s][LocationType::Left]  = new double[ny];
    buffer_map[s][LocationType::Right] = new double[ny];
    buffer_map[s][LocationType::Down]  = new double[snx];
    buffer_map[s][LocationType::Up]    = new double[snx];

    position_type = VariablePositionType::Center;
}

void Variable2DSlabX::set_boundary_value(Domain2DUniformMPI* s, LocationType loc, double in_value)
{
    check_geometry(s);

    int level = 0;
    if (slab_parent_to_level.find(s->get_uuid()) != slab_parent_to_level.end())
        level = slab_parent_to_level[s->get_uuid()];
    else
        return;

    int snx = hierarchical_slab_nxs[level];
    int ny  = s->ny;

    has_boundary_value_map[s][loc] = true;
    if (loc == LocationType::Left || loc == LocationType::Right)
    {
        boundary_value_map[s][loc] = new double[ny];
        for (int j = 0; j < ny; j++)
            boundary_value_map[s][loc][j] = in_value;

        if (loc == LocationType::Left)
            left_up_corner_value_map[s] = in_value;
    }
    else if (loc == LocationType::Down || loc == LocationType::Up)
    {
        boundary_value_map[s][loc] = new double[snx];
        for (int i = 0; i < snx; i++)
            boundary_value_map[s][loc][i] = in_value;

        if (loc == LocationType::Down)
            right_down_corner_value_map[s] = in_value;
    }
}

void Variable2DSlabX::set_boundary_value_from_func_global(Domain2DUniformMPI*                   s,
                                                          LocationType                          loc,
                                                          std::function<double(double, double)> f)
{
    check_geometry(s);

    int level = 0;
    if (slab_parent_to_level.find(s->get_uuid()) != slab_parent_to_level.end())
        level = slab_parent_to_level[s->get_uuid()];
    else
        return;

    int slab_rank = hierarchical_slab_ranks[level];
    int slab_size = hierarchical_slab_sizes[level];

    // Left bound at rank=0, Right bound at rank=size-1
    if (loc == LocationType::Left || loc == LocationType::Right)
        if (slab_rank != 0 && slab_rank != slab_size - 1)
            return;

    int nx_slab = hierarchical_slab_nxs[level];
    int nx_disp = hierarchical_slab_disps[level];

    has_boundary_value_map[s][loc] = true;

    double shift_x = 0.5;
    double shift_y = 0.5;

    switch (position_type)
    {
        case VariablePositionType::Center:
            shift_x = 0.5;
            shift_y = 0.5;
            break;
        case VariablePositionType::XFace:
            shift_x = 0.0;
            shift_y = 0.5;
            break;
        case VariablePositionType::YFace:
            shift_x = 0.5;
            shift_y = 0.0;
            break;
        case VariablePositionType::Corner:
            shift_x = 0.0;
            shift_y = 0.0;
            break;
        default:
            break;
    }

    if (loc == LocationType::Left || loc == LocationType::Right)
    {
        int j_size;
        if (position_type == VariablePositionType::Corner)
            j_size = s->ny + 1;
        else
            j_size = s->ny;

        boundary_value_map[s][loc] = new double[j_size];

        // Determine i index (ghost node or boundary node)
        // Left: i = -1. Right: i = nx.
        int i_idx = (loc == LocationType::Left) ? -1 : s->nx;

        for (int j = 0; j < j_size; j++)
        {
            double gx                     = s->get_offset_x() + (shift_x + i_idx) * s->get_hx();
            double gy                     = s->get_offset_y() + (shift_y + j) * s->get_hy();
            boundary_value_map[s][loc][j] = f(gx, gy);
        }
    }
    else if (loc == LocationType::Down || loc == LocationType::Up)
    {
        int i_size = nx_slab;

        boundary_value_map[s][loc] = new double[i_size];

        // Down: j = -1. Up: j = ny.
        int j_idx = (loc == LocationType::Down) ? -1 : s->ny;

        for (int i = 0; i < i_size; i++)
        {
            double gx                     = s->get_offset_x() + (shift_x + (i + nx_disp)) * s->get_hx();
            double gy                     = s->get_offset_y() + (shift_y + j_idx) * s->get_hy();
            boundary_value_map[s][loc][i] = f(gx, gy);
        }
    }
}

void Variable2DSlabX::fill_boundary_value_from_func_global(std::function<double(double, double)> f)
{
    if (geometry == nullptr)
        throw std::runtime_error("Variable2DSlabX has no geometry set");

    for (auto kv : field_map)
    {
        Domain2DUniform* domain = kv.first;

        auto& type_map = boundary_type_map[domain];
        auto& has_map  = has_boundary_value_map[domain];
        for (LocationType loc : kBoundaryLocations2D)
        {
            auto type_it = type_map.find(loc);
            if (type_it == type_map.end())
                continue;
            if (type_it->second == PDEBoundaryType::Adjacented || type_it->second == PDEBoundaryType::Null)
                continue;
            auto has_it = has_map.find(loc);
            if (has_it != has_map.end() && has_it->second)
                continue;

            set_boundary_value_from_func_global(static_cast<Domain2DUniformMPI*>(domain), loc, f);
        }
    }
}

void Variable2DSlabX::set_value_from_func_global(std::function<double(double, double)> func)
{
    double shift_x = 0.5;
    double shift_y = 0.5;

    switch (position_type)
    {
        case VariablePositionType::Center:
            shift_x = 0.5;
            shift_y = 0.5;
            break;
        case VariablePositionType::XFace:
            shift_x = 0.0;
            shift_y = 0.5;
            break;
        case VariablePositionType::YFace:
            shift_x = 0.5;
            shift_y = 0.0;
            break;
        case VariablePositionType::Corner:
            shift_x = 0.0;
            shift_y = 0.0;
            break;
        default:
            // Or throw exception? Defaulting to Center for now or Null
            break;
    }

    for (auto& pair : field_map)
    {
        Domain2DUniformMPI* s = static_cast<Domain2DUniformMPI*>(pair.first);
        field2*             f = pair.second;

        int level   = slab_parent_to_level[s->get_uuid()];
        int nx_disp = hierarchical_slab_disps[level];

        OPENMP_PARALLEL_FOR()
        for (int i = 0; i < f->get_nx(); ++i)
        {
            for (int j = 0; j < f->get_ny(); ++j)
            {
                double gx  = s->get_offset_x() + (shift_x + (i + nx_disp)) * s->get_hx();
                double gy  = s->get_offset_y() + (shift_y + j) * s->get_hy();
                (*f)(i, j) = func(gx, gy);
            }
        }
    }
}

void Variable2DSlabX::print_slab_info()
{
    if (mpi_rank == 0)
    {
        std::cout << "Variable2DSlabX slab info:" << std::endl;
        std::cout << "PositionType = " << position_type << std::endl;
        int level_num = 0;
        for (auto level : geometry->hierarchical_solve_levels)
        {
            std::cout << "Level " << level_num << std::endl;
            level_num++;

            for (int i = 0; i < level.size(); i++)
            {
                Domain2DUniformMPI* domain = static_cast<Domain2DUniformMPI*>(level[i]);
                std::cout << " domain name = " << domain->name << " uuid = " << domain->get_uuid();
                std::cout << " [" << domain->get_nx() << ", " << domain->get_ny() << "]" << std::endl;
            }
        }
    }

    for (int level = 0; level < hierarchical_slab_parents.size(); level++)
    {
        if (mpi_rank == 0)
            std::cout << "Level " << level << std::endl;

        for (int i = 0; i < mpi_size; i++)
        {
            if (i == mpi_rank)
            {
                Domain2DUniformMPI* domain = static_cast<Domain2DUniformMPI*>(hierarchical_slab_parents[level]);
                std::cout << "rank " << mpi_rank;
                std::cout << " domain " << domain->get_uuid();
                std::cout << " [" << hierarchical_slab_nxs[level] << ", " << field_map[domain]->get_ny() << "] + "
                          << hierarchical_slab_disps[level] << std::endl;
            }
            MPI_Barrier(MPI_COMM_WORLD);
        }
    }
}

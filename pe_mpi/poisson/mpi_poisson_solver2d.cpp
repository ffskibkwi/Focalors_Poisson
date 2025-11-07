#include "mpi_poisson_solver2d.h"

#include <cstring>

MPIPoissonSolver2D::MPIPoissonSolver2D(int                in_nx,
                                       int                in_ny,
                                       double             in_hx,
                                       double             in_hy,
                                       PDEBoundaryType    in_boundary_type_left,
                                       PDEBoundaryType    in_boundary_type_right,
                                       PDEBoundaryType    in_boundary_type_down,
                                       PDEBoundaryType    in_boundary_type_up,
                                       int                in_num_proc,
                                       EnvironmentConfig* in_env_config,
                                       MPI_Comm           in_comm)
    : nx(in_nx)
    , ny(in_ny)
    , hx(in_hx)
    , hy(in_hy)
    , boundary_type_left(in_boundary_type_left)
    , boundary_type_right(in_boundary_type_right)
    , boundary_type_down(in_boundary_type_down)
    , boundary_type_up(in_boundary_type_up)
    , num_proc(in_num_proc)
{
    env_config = in_env_config;
    setup_comm(in_comm);
    setup_problem_params();
    build_decomposition();
    build_local_solvers();
}

MPIPoissonSolver2D::MPIPoissonSolver2D(Domain2DUniform* in_domain,
                                       Variable*        in_variable,
                                       int              in_num_proc,
                                       EnvironmentConfig* in_env_config,
                                       MPI_Comm           in_comm)
    : domain(in_domain)
    , var(in_variable)
    , nx(in_domain->nx)
    , ny(in_domain->ny)
    , hx(in_domain->hx)
    , hy(in_domain->hy)
    , num_proc(in_num_proc)
{
    env_config           = in_env_config;
    boundary_type_left   = var->boundary_type_map[domain][LocationType::Left];
    boundary_type_right  = var->boundary_type_map[domain][LocationType::Right];
    boundary_type_down   = var->boundary_type_map[domain][LocationType::Down];
    boundary_type_up     = var->boundary_type_map[domain][LocationType::Up];

    setup_comm(in_comm);
    setup_problem_params();
    build_decomposition();
    build_local_solvers();
}

MPIPoissonSolver2D::~MPIPoissonSolver2D()
{
    release_local_solvers();
    if (active_comm != MPI_COMM_NULL)
    {
        MPI_Comm_free(&active_comm);
        active_comm = MPI_COMM_NULL;
    }
}

void MPIPoissonSolver2D::setup_comm(MPI_Comm in_comm)
{
    world_comm = in_comm;
    MPI_Comm_rank(world_comm, &world_rank);
    MPI_Comm_size(world_comm, &world_size);

    int desired = num_proc;
    if (desired <= 0 || desired > world_size)
        desired = world_size;

    int color = (world_rank < desired) ? 1 : MPI_UNDEFINED;
    MPI_Comm_split(world_comm, color, world_rank, &active_comm);

    if (color == 1)
    {
        is_active = true;
        MPI_Comm_rank(active_comm, &active_rank);
        MPI_Comm_size(active_comm, &active_size);
    }
    else
    {
        is_active = false;
        active_rank = -1;
        active_size = 0;
    }
}

void MPIPoissonSolver2D::setup_problem_params()
{
    // nothing else for now
}

void MPIPoissonSolver2D::build_decomposition()
{
    if (!is_active)
        return;

    split_1d(nx, active_size, i_counts, i_displs);
    split_1d(ny, active_size, j_counts, j_displs);

    local_i_count = i_counts[active_rank];
    local_j_count = j_counts[active_rank];
}

void MPIPoissonSolver2D::build_local_solvers()
{
    if (!is_active)
        return;

    poisson_fft_y = new PoissonFFT2D_PP();
    auto isDirLike = [](PDEBoundaryType t) {
        return t == PDEBoundaryType::Dirichlet || t == PDEBoundaryType::Adjacented;
    };

    if (boundary_type_down == PDEBoundaryType::Periodic && boundary_type_up == PDEBoundaryType::Periodic)
    {
        poisson_fft_y = new PoissonFFT2D_PP();
    }
    else if (boundary_type_down == PDEBoundaryType::Neumann && boundary_type_up == PDEBoundaryType::Neumann)
    {
        poisson_fft_y = new PoissonFFT2D_NN();
    }
    else if (isDirLike(boundary_type_down) && isDirLike(boundary_type_up))
    {
        poisson_fft_y = new PoissonFFT2D_DD();
    }
    else if ((isDirLike(boundary_type_down) && boundary_type_up == PDEBoundaryType::Neumann) ||
             (boundary_type_down == PDEBoundaryType::Neumann && isDirLike(boundary_type_up)))
    {
        poisson_fft_y = new PoissonFFT2D_DN();
    }
    poisson_fft_y->init(local_i_count, ny);

    build_local_x_diag();

    bool is_no_Dirichlet = !(isDirLike(boundary_type_left) || isDirLike(boundary_type_right) ||
                             isDirLike(boundary_type_down) || isDirLike(boundary_type_up));
    bool has_last_vector = true;

    chasing_method_x = new ChasingMethod2D();
    // For transposed layout: (local_j_count, nx)
    chasing_method_x->init(local_j_count, nx, local_x_diag, is_no_Dirichlet, has_last_vector,
                           boundary_type_left, boundary_type_right);
}

void MPIPoissonSolver2D::release_local_solvers()
{
    if (local_x_diag)
    {
        delete[] local_x_diag;
        local_x_diag = nullptr;
    }
    if (poisson_fft_y)
    {
        delete poisson_fft_y;
        poisson_fft_y = nullptr;
    }
    if (chasing_method_x)
    {
        delete chasing_method_x;
        chasing_method_x = nullptr;
    }
}

void MPIPoissonSolver2D::cal_lambda(std::vector<double>& lambda_y_out) const
{
    lambda_y_out.resize(ny);
    for (int i = 1; i <= ny; i++)
    {
        if (boundary_type_down == PDEBoundaryType::Periodic && boundary_type_up == PDEBoundaryType::Periodic)
        {
            lambda_y_out[i - 1] = -2.0 + 2.0 * std::cos(2.0 * pi / ny * std::floor(i / 2.0));
        }
        else if (boundary_type_down == PDEBoundaryType::Neumann && boundary_type_up == PDEBoundaryType::Neumann)
        {
            lambda_y_out[i - 1] = -2.0 + 2.0 * std::cos(pi / ny * i);
        }
        else if ((boundary_type_down == PDEBoundaryType::Dirichlet || boundary_type_down == PDEBoundaryType::Adjacented) &&
                 (boundary_type_up == PDEBoundaryType::Dirichlet || boundary_type_up == PDEBoundaryType::Adjacented))
        {
            lambda_y_out[i - 1] = -2.0 + 2.0 * std::cos(pi / (ny + 1) * i);
        }
        else if (((boundary_type_down == PDEBoundaryType::Dirichlet || boundary_type_down == PDEBoundaryType::Adjacented) &&
                  boundary_type_up == PDEBoundaryType::Neumann) ||
                 (boundary_type_down == PDEBoundaryType::Neumann &&
                  (boundary_type_up == PDEBoundaryType::Dirichlet || boundary_type_up == PDEBoundaryType::Adjacented)))
        {
            lambda_y_out[i - 1] = -2.0 + 2.0 * std::cos(2.0 * pi * i / (2 * ny + 1));
        }
    }
}

void MPIPoissonSolver2D::build_local_x_diag()
{
    if (!is_active)
        return;

    std::vector<double> lambda_y_vec;
    cal_lambda(lambda_y_vec);

    // Build local slice for owned j indices
    local_x_diag = new double[local_j_count];
    int j_start  = j_displs[active_rank];
    for (int j = 0; j < local_j_count; ++j)
    {
        local_x_diag[j] = -2.0 + hx * hx / hy / hy * lambda_y_vec[j_start + j];
    }
}

void MPIPoissonSolver2D::boundary_assembly(field2& f)
{
    auto& var_has_map   = var->has_boundary_value_map[domain];
    auto& var_value_map = var->boundary_value_map[domain];

    if (boundary_type_left == PDEBoundaryType::Dirichlet && var_has_map[LocationType::Left])
    {
        double* boundary_value = var_value_map[LocationType::Left];
        for (int j = 0; j < ny; j++)
            f(0, j) -= boundary_value[j] / hx / hx;
    }
    if (boundary_type_right == PDEBoundaryType::Dirichlet && var_has_map[LocationType::Right])
    {
        double* boundary_value = var_value_map[LocationType::Right];
        for (int j = 0; j < ny; j++)
            f(nx - 1, j) -= boundary_value[j] / hx / hx;
    }

    if (boundary_type_down == PDEBoundaryType::Dirichlet && var_has_map[LocationType::Down])
    {
        double* boundary_value = var_value_map[LocationType::Down];
        for (int i = 0; i < nx; i++)
            f(i, 0) -= boundary_value[i] / hy / hy;
    }
    if (boundary_type_up == PDEBoundaryType::Dirichlet && var_has_map[LocationType::Up])
    {
        double* boundary_value = var_value_map[LocationType::Up];
        for (int i = 0; i < nx; i++)
            f(i, ny - 1) -= boundary_value[i] / hy / hy;
    }

    if (boundary_type_left == PDEBoundaryType::Neumann && var_has_map[LocationType::Left])
    {
        double* boundary_value = var_value_map[LocationType::Left];
        for (int j = 0; j < ny; j++)
            f(0, j) += boundary_value[j] / hx;
    }
    if (boundary_type_right == PDEBoundaryType::Neumann && var_has_map[LocationType::Right])
    {
        double* boundary_value = var_value_map[LocationType::Left];
        for (int j = 0; j < ny; j++)
            f(nx - 1, j) -= boundary_value[j] / hx;
    }
    if (boundary_type_down == PDEBoundaryType::Neumann && var_has_map[LocationType::Down])
    {
        double* boundary_value = var_value_map[LocationType::Down];
        for (int i = 0; i < nx; i++)
            f(i, 0) += boundary_value[i] / hy;
    }
    if (boundary_type_up == PDEBoundaryType::Neumann && var_has_map[LocationType::Up])
    {
        double* boundary_value = var_value_map[LocationType::Up];
        for (int i = 0; i < nx; i++)
            f(i, ny - 1) -= boundary_value[i] / hy;
    }
}

void MPIPoissonSolver2D::scatter_global_f_to_local(const field2& f_global, field2& f_local) const
{
    // Allocate local
    f_local.init(local_i_count, ny, "local_f");

    // Prepare counts/displs in elements
    std::vector<int> sendcounts(active_size), displs(active_size);
    for (int r = 0; r < active_size; ++r)
    {
        sendcounts[r] = i_counts[r] * ny;
        displs[r]     = i_displs[r] * ny;
    }

    MPI_Scatterv(
        f_global.value,
        sendcounts.data(),
        displs.data(),
        MPI_DOUBLE,
        f_local.value,
        sendcounts[active_rank],
        MPI_DOUBLE,
        0,
        active_comm);
}

void MPIPoissonSolver2D::gather_local_to_global_f(const field2& f_local, field2& f_global) const
{
    std::vector<int> recvcounts(active_size), displs(active_size);
    for (int r = 0; r < active_size; ++r)
    {
        recvcounts[r] = i_counts[r] * ny;
        displs[r]     = i_displs[r] * ny;
    }

    MPI_Gatherv(
        f_local.value,
        recvcounts[active_rank],
        MPI_DOUBLE,
        f_global.value,
        recvcounts.data(),
        displs.data(),
        MPI_DOUBLE,
        0,
        active_comm);
}

void MPIPoissonSolver2D::distributed_transpose_i_to_j(const field2& in_i_slab, field2& out_j_slab) const
{
    // in_i_slab: (local_i_count, ny)
    // out_j_slab: (local_j_count, nx)

    // Pack send buffers per-destination rank
    std::vector<int> sendcounts(active_size), sdispls(active_size);
    std::vector<int> recvcounts(active_size), rdispls(active_size);
    for (int r = 0; r < active_size; ++r)
    {
        sendcounts[r] = local_i_count * j_counts[r];
        recvcounts[r] = i_counts[r] * local_j_count;
    }
    sdispls[0] = 0;
    rdispls[0] = 0;
    for (int r = 1; r < active_size; ++r)
    {
        sdispls[r] = sdispls[r - 1] + sendcounts[r - 1];
        rdispls[r] = rdispls[r - 1] + recvcounts[r - 1];
    }

    std::vector<double> sendbuf(sdispls.back() + sendcounts.back());
    std::vector<double> recvbuf(rdispls.back() + recvcounts.back());

    for (int r = 0; r < active_size; ++r)
    {
        int j0 = j_displs[r];
        int jc = j_counts[r];
        double* dst = sendbuf.data() + sdispls[r];
        for (int i = 0; i < local_i_count; ++i)
        {
            std::memcpy(dst + i * jc, in_i_slab.get_ptr(i, j0), sizeof(double) * jc);
        }
    }

    MPI_Alltoallv(sendbuf.data(), sendcounts.data(), sdispls.data(), MPI_DOUBLE,
                  recvbuf.data(), recvcounts.data(), rdispls.data(), MPI_DOUBLE,
                  active_comm);

    // Unpack into out_j_slab (local_j_count, nx)
    out_j_slab.init(local_j_count, nx, "local_transposed");
    for (int s = 0; s < active_size; ++s)
    {
        int ic = i_counts[s];
        int i0 = i_displs[s];
        const double* src = recvbuf.data() + rdispls[s];
        for (int il = 0; il < ic; ++il)
        {
            int i_global = i0 + il;
            for (int jl = 0; jl < local_j_count; ++jl)
            {
                double v = src[il * local_j_count + jl];
                out_j_slab(jl, i_global) = v;
            }
        }
    }
}

void MPIPoissonSolver2D::distributed_transpose_j_to_i(const field2& in_j_slab, field2& out_i_slab) const
{
    // in_j_slab: (local_j_count, nx)
    // out_i_slab: (local_i_count, ny)

    std::vector<int> sendcounts(active_size), sdispls(active_size);
    std::vector<int> recvcounts(active_size), rdispls(active_size);
    for (int r = 0; r < active_size; ++r)
    {
        sendcounts[r] = i_counts[r] * local_j_count;
        recvcounts[r] = local_i_count * j_counts[r];
    }
    sdispls[0] = 0;
    rdispls[0] = 0;
    for (int r = 1; r < active_size; ++r)
    {
        sdispls[r] = sdispls[r - 1] + sendcounts[r - 1];
        rdispls[r] = rdispls[r - 1] + recvcounts[r - 1];
    }

    std::vector<double> sendbuf(sdispls.back() + sendcounts.back());
    std::vector<double> recvbuf(rdispls.back() + recvcounts.back());

    for (int r = 0; r < active_size; ++r)
    {
        int i0 = i_displs[r];
        int ic = i_counts[r];
        double* dst = sendbuf.data() + sdispls[r];
        for (int il = 0; il < ic; ++il)
        {
            int i_global = i0 + il;
            for (int jl = 0; jl < local_j_count; ++jl)
            {
                dst[il * local_j_count + jl] = in_j_slab(jl, i_global);
            }
        }
    }

    MPI_Alltoallv(sendbuf.data(), sendcounts.data(), sdispls.data(), MPI_DOUBLE,
                  recvbuf.data(), recvcounts.data(), rdispls.data(), MPI_DOUBLE,
                  active_comm);

    out_i_slab.init(local_i_count, ny, "local_back");
    for (int s = 0; s < active_size; ++s)
    {
        int jc = j_counts[s];
        int j0 = j_displs[s];
        const double* src = recvbuf.data() + rdispls[s];
        for (int il = 0; il < local_i_count; ++il)
        {
            for (int jl = 0; jl < jc; ++jl)
            {
                out_i_slab(il, j0 + jl) = src[il * jc + jl];
            }
        }
    }
}

void MPIPoissonSolver2D::solve(field2& f)
{
    // Only active ranks participate
    if (!is_active)
        return;

    if (env_config && env_config->showCurrentStep && active_rank == 0)
        std::cout << "[Poisson][MPI] solve: start" << std::endl;

    // 1) Boundary assembly on root, since full field is only on rank 0
    if (active_rank == 0 && var && domain)
    {
        boundary_assembly(f);
    }

    // 2) Scatter to i-slab local fields
    field2 f_local;
    scatter_global_f_to_local(f, f_local);

    // 3) Y-direction transform locally: f_local -> fhat_local
    field2 fhat_local;
    fhat_local.init(local_i_count, ny, "fhat_local");
    poisson_fft_y->transform(f_local, fhat_local);

    // 4) Global transpose: i-slab -> j-slab
    field2 fhat_local_T;
    distributed_transpose_i_to_j(fhat_local, fhat_local_T);

    // 5) X-direction chasing on transposed slabs (local_j_count, nx)
    field2 phat_local_T;
    phat_local_T.init(local_j_count, nx, "phat_local_T");
    chasing_method_x->chasing(fhat_local_T, phat_local_T);

    // 6) Global transpose back: j-slab -> i-slab
    field2 phat_local;
    distributed_transpose_j_to_i(phat_local_T, phat_local);

    // 7) Inverse Y-transform locally
    field2 p_local;
    p_local.init(local_i_count, ny, "p_local");
    poisson_fft_y->transform_transpose(phat_local, p_local);

    // 8) Gather back to root
    gather_local_to_global_f(p_local, f);

    if (env_config && env_config->showCurrentStep && active_rank == 0)
        std::cout << "[Poisson][MPI] solve: done" << std::endl;
}

void MPIPoissonSolver2D::split_1d(int n, int p, std::vector<int>& counts, std::vector<int>& displs)
{
    counts.assign(p, 0);
    displs.assign(p, 0);
    int base = n / p;
    int rem  = n % p;
    int off  = 0;
    for (int r = 0; r < p; ++r)
    {
        counts[r] = base + (r < rem ? 1 : 0);
        displs[r] = off;
        off += counts[r];
    }
}



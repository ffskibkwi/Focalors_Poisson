#include "poisson_solver2d_slab_x.h"

#include "base/parallel/mpi/mpi_misc.h"
#include "base/parallel/mpi/transpose_slab.h"

PoissonSolver2DSlabX::PoissonSolver2DSlabX(int             in_nx,
                                           int             in_ny,
                                           double          in_hx,
                                           double          in_hy,
                                           PDEBoundaryType in_boundary_type_left,
                                           PDEBoundaryType in_boundary_type_right,
                                           PDEBoundaryType in_boundary_type_down,
                                           PDEBoundaryType in_boundary_type_up,
                                           MPI_Comm        in_communicator)
    : PoissonSolver2DBase(in_nx,
                          in_ny,
                          in_hx,
                          in_hy,
                          in_boundary_type_left,
                          in_boundary_type_right,
                          in_boundary_type_down,
                          in_boundary_type_up)
    , communicator(in_communicator)
{
    init();
}

PoissonSolver2DSlabX::PoissonSolver2DSlabX(Domain2DUniform* in_domain,
                                           Variable2D*      in_variable,
                                           MPI_Comm         in_communicator)
    : PoissonSolver2DBase(in_domain->nx,
                          in_domain->ny,
                          in_domain->hx,
                          in_domain->hy,
                          in_variable->boundary_type_map[in_domain][LocationType::Left],
                          in_variable->boundary_type_map[in_domain][LocationType::Right],
                          in_variable->boundary_type_map[in_domain][LocationType::Down],
                          in_variable->boundary_type_map[in_domain][LocationType::Up])
    , domain_name(in_domain->name)
    , communicator(in_communicator)
{
    init();
}

void PoissonSolver2DSlabX::init()
{
    MPI_Comm_rank(communicator, &mpi_rank);
    MPI_Comm_size(communicator, &mpi_size);

    snx     = MPIUtils::get_slab_length(nx, mpi_rank, mpi_size);
    sny     = MPIUtils::get_slab_length(ny, mpi_rank, mpi_size);
    ny_disp = MPIUtils::get_slab_displacement(ny, mpi_rank, mpi_size);

    f_hat.init(snx, ny, "f_hat");
    f_hat_T.init(sny, nx, "f_hat_T");
    p_hat_T.init(sny, nx, "p_hat_T");
    p_hat.init(snx, ny, "p_hat");

    create_fft(poisson_fft_y, boundary_type_down, boundary_type_up, snx, ny);

    double* lambda_y = new double[sny];
    x_diag           = new double[sny];

    cal_lambda(lambda_y, ny, 1 + ny_disp, sny, boundary_type_down, boundary_type_up);

    for (int j = 0; j < sny; j++)
    {
        x_diag[j] = -2.0 + hx * hx / hy / hy * (lambda_y[j] - 2.0);
    }
    delete[] lambda_y;

    bool is_no_Dirichlet = !(isDirLike(boundary_type_left) || isDirLike(boundary_type_right) ||
                             isDirLike(boundary_type_down) || isDirLike(boundary_type_up));
    bool has_last_vector = mpi_rank == mpi_size - 1;

    chasing_method_x = new ChasingMethod2D();
    chasing_method_x->init(sny, nx, x_diag, is_no_Dirichlet, has_last_vector, boundary_type_left, boundary_type_right);

    buf_slab_x = new double[snx * ny];
    buf_slab_y = new double[sny * nx];
}

PoissonSolver2DSlabX::~PoissonSolver2DSlabX()
{
    delete[] buf_slab_y;
    delete[] buf_slab_x;
}

void PoissonSolver2DSlabX::solve(field2& f)
{
    EnvironmentConfig& env_cfg = EnvironmentConfig::Get();

    SCOPE_TIMER("PoissonSolver2DSlabX::solve", TimeRecordType::None, env_cfg.track_pe_solve_detail_time);

    if (env_cfg.showCurrentStep)
        std::cout << "[Poisson] solve: start" << std::endl;

    poisson_fft_y->transform(f, f_hat);

    MPIUtils::transpose_2d_slab_sync(f_hat, buf_slab_x, f_hat_T, buf_slab_y, communicator);

    chasing_method_x->chasing(f_hat_T, p_hat_T);

    MPIUtils::transpose_2d_slab_sync(p_hat_T, buf_slab_y, p_hat, buf_slab_x, communicator);

    poisson_fft_y->transform_transpose(p_hat, f);

    solve_call_count++;
    if (env_cfg.showCurrentStep)
    {
        double s_f = f.sum();
        std::cout << "[Poisson] solve: done, f.sum=" << s_f << std::endl;
    }
}

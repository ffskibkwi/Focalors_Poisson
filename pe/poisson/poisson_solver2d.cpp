#include "poisson_solver2d.h"

PoissonSolver2D::PoissonSolver2D(int             in_nx,
                                 int             in_ny,
                                 double          in_hx,
                                 double          in_hy,
                                 PDEBoundaryType in_boundary_type_left,
                                 PDEBoundaryType in_boundary_type_right,
                                 PDEBoundaryType in_boundary_type_down,
                                 PDEBoundaryType in_boundary_type_up)
    : PoissonSolver2DBase(in_nx,
                          in_ny,
                          in_hx,
                          in_hy,
                          in_boundary_type_left,
                          in_boundary_type_right,
                          in_boundary_type_down,
                          in_boundary_type_up)
{
    init();
}

PoissonSolver2D::PoissonSolver2D(Domain2DUniform* in_domain, Variable2D* in_variable)
    : PoissonSolver2DBase(in_domain->nx,
                          in_domain->ny,
                          in_domain->hx,
                          in_domain->hy,
                          in_variable->boundary_type_map[in_domain][LocationType::Left],
                          in_variable->boundary_type_map[in_domain][LocationType::Right],
                          in_variable->boundary_type_map[in_domain][LocationType::Down],
                          in_variable->boundary_type_map[in_domain][LocationType::Up])
    , domain_name(in_domain->name)
{
    init();
}

void PoissonSolver2D::init()
{
    buffer.init(nx, ny, "buffer");

    create_fft(poisson_fft_y, boundary_type_down, boundary_type_up, nx, ny);

    double* lambda_y = new double[ny];
    x_diag           = new double[ny];

    cal_lambda(lambda_y, ny, 1, ny, boundary_type_down, boundary_type_up);

    for (int j = 0; j < ny; j++)
    {
        x_diag[j] = -2.0 + hx * hx / hy / hy * (lambda_y[j] - 2.0);
    }
    delete[] lambda_y;

    bool is_no_Dirichlet = !(isDirLike(boundary_type_left) || isDirLike(boundary_type_right) ||
                             isDirLike(boundary_type_down) || isDirLike(boundary_type_up));
    bool has_last_vector = true;

    chasing_method_x = new ChasingMethod2D();
    chasing_method_x->init(ny, nx, x_diag, is_no_Dirichlet, has_last_vector, boundary_type_left, boundary_type_right);
}

void PoissonSolver2D::solve(field2& f)
{
    EnvironmentConfig& env_cfg = EnvironmentConfig::Get();

    SCOPE_TIMER("PoissonSolver2D::solve", TimeRecordType::None, env_cfg.track_pe_solve_detail_time);

    if (env_cfg.showCurrentStep)
        std::cout << "[Poisson] solve: start" << std::endl;

    buffer.set_size(nx, ny);
    poisson_fft_y->transform(f, buffer);

    buffer.transpose(f);
    buffer.set_size(ny, nx);
    chasing_method_x->chasing(f, buffer);

    buffer.transpose(f);
    buffer.set_size(nx, ny);
    poisson_fft_y->transform_transpose(f, buffer);

    std::swap(f, buffer);

    solve_call_count++;
    if (env_cfg.showCurrentStep)
    {
        double s_f = f.sum();
        std::cout << "[Poisson] solve: done, f.sum=" << s_f << std::endl;
    }
}

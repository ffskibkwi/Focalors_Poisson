#include "poisson_solver2d.h"

PoissonSolver2D::PoissonSolver2D(int             in_nx,
                                 int             in_ny,
                                 double          in_hx,
                                 double          in_hy,
                                 PDEBoundaryType in_boundary_type_xneg,
                                 PDEBoundaryType in_boundary_type_xpos,
                                 PDEBoundaryType in_boundary_type_yneg,
                                 PDEBoundaryType in_boundary_type_ypos)
    : PoissonSolver2DBase(in_nx,
                          in_ny,
                          in_hx,
                          in_hy,
                          in_boundary_type_xneg,
                          in_boundary_type_xpos,
                          in_boundary_type_yneg,
                          in_boundary_type_ypos)
{
    init();
}

PoissonSolver2D::PoissonSolver2D(Domain2DUniform* in_domain, Variable2D* in_variable)
    : PoissonSolver2DBase(in_domain->nx,
                          in_domain->ny,
                          in_domain->hx,
                          in_domain->hy,
                          in_variable->boundary_type_map[in_domain][LocationType::XNegative],
                          in_variable->boundary_type_map[in_domain][LocationType::XPositive],
                          in_variable->boundary_type_map[in_domain][LocationType::YNegative],
                          in_variable->boundary_type_map[in_domain][LocationType::YPositive])
    , domain_name(in_domain->name)
{
    init();
}

void PoissonSolver2D::init()
{
    buffer.init(nx, ny, "buffer");

    create_fft(poisson_fft_y, boundary_type_yneg, boundary_type_ypos, nx, ny);

    double* lambda_y = new double[ny];
    x_diag           = new double[ny];

    cal_lambda(lambda_y, ny, 1, ny, boundary_type_yneg, boundary_type_ypos);

    for (int j = 0; j < ny; j++)
    {
        x_diag[j] = -2.0 + hx * hx / hy / hy * (lambda_y[j] - 2.0);
    }
    delete[] lambda_y;

    bool is_no_Dirichlet = !(isDirLike(boundary_type_xneg) || isDirLike(boundary_type_xpos) ||
                             isDirLike(boundary_type_yneg) || isDirLike(boundary_type_ypos));
    bool has_last_vector = true;

    chasing_method_x = new ChasingMethod2D();
    chasing_method_x->init(ny, nx, x_diag, is_no_Dirichlet, has_last_vector, boundary_type_xneg, boundary_type_xpos);
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

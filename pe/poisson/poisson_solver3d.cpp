#include "poisson_solver3d.h"

#include "instrumentor/timer.h"

PoissonSolver3D::PoissonSolver3D(int             in_nx,
                                 int             in_ny,
                                 int             in_nz,
                                 double          in_hx,
                                 double          in_hy,
                                 double          in_hz,
                                 PDEBoundaryType in_boundary_type_left,
                                 PDEBoundaryType in_boundary_type_right,
                                 PDEBoundaryType in_boundary_type_front,
                                 PDEBoundaryType in_boundary_type_back,
                                 PDEBoundaryType in_boundary_type_down,
                                 PDEBoundaryType in_boundary_type_up)
    : PoissonSolver3DBase(in_nx,
                          in_ny,
                          in_nz,
                          in_hx,
                          in_hy,
                          in_hz,
                          in_boundary_type_left,
                          in_boundary_type_right,
                          in_boundary_type_front,
                          in_boundary_type_back,
                          in_boundary_type_down,
                          in_boundary_type_up)
{
    init();
}

PoissonSolver3D::PoissonSolver3D(Domain3DUniform* in_domain, Variable3D* in_variable)
    : PoissonSolver3DBase(in_domain->nx,
                          in_domain->ny,
                          in_domain->nz,
                          in_domain->hx,
                          in_domain->hy,
                          in_domain->hz,
                          in_variable->boundary_type_map[in_domain][LocationType::Left],
                          in_variable->boundary_type_map[in_domain][LocationType::Right],
                          in_variable->boundary_type_map[in_domain][LocationType::Front],
                          in_variable->boundary_type_map[in_domain][LocationType::Back],
                          in_variable->boundary_type_map[in_domain][LocationType::Down],
                          in_variable->boundary_type_map[in_domain][LocationType::Up])
    , domain_name(in_domain->name)
{
    init();
}

void PoissonSolver3D::init()
{
    buffer.init(nx, ny, nz, "buffer");

    create_fft(poisson_fft_z, boundary_type_down, boundary_type_up, nx, ny, nz);
    create_fft(poisson_fft_y, boundary_type_front, boundary_type_back, nz, nx, ny);

    double* lambda_z = new double[nz];
    double* lambda_y = new double[ny];
    x_diag           = new double*[nz];

    cal_lambda(lambda_z, nz, 1, nz, boundary_type_down, boundary_type_up);
    cal_lambda(lambda_y, ny, 1, ny, boundary_type_front, boundary_type_back);

    for (int k = 0; k < nz; k++)
    {
        x_diag[k] = new double[ny];
        for (int j = 0; j < ny; j++)
        {
            x_diag[k][j] = -2.0 + hx * hx / hy / hy * lambda_y[j] + hx * hx / hz / hz * lambda_z[k];
        }
    }
    delete[] lambda_y;
    delete[] lambda_z;

    bool is_no_Dirichlet =
        !(isDirLike(boundary_type_left) || isDirLike(boundary_type_right) || isDirLike(boundary_type_front) ||
          isDirLike(boundary_type_back) || isDirLike(boundary_type_down) || isDirLike(boundary_type_up));
    bool has_last_vector = true;

    chasing_method_x = new ChasingMethod3D();
    chasing_method_x->init(
        nz, ny, nx, x_diag, is_no_Dirichlet, has_last_vector, boundary_type_left, boundary_type_right);
}

void PoissonSolver3D::solve(field3& f)
{
    EnvironmentConfig& env_cfg = EnvironmentConfig::Get();

    SCOPE_TIMER("PoissonSolver3D::solve", TimeRecordType::None, env_cfg.track_pe_solve_detail_time);
    SCOPE_TIMER(env_cfg.pe_solve_total_name, TimeRecordType::Accumulate, false);

    if (env_cfg.showCurrentStep)
        std::cout << "[Poisson] solve: start (domain " << domain_name << ")" << std::endl;

    buffer.set_size(nx, ny, nz);
    poisson_fft_z->transform(f, buffer);

    f.set_size(nz, nx, ny);
    buffer.transpose(f, {2, 0, 1});

    buffer.set_size(nz, nx, ny);
    poisson_fft_y->transform(f, buffer);

    f.set_size(ny, nz, nx);
    buffer.transpose(f, {0, 2, 1});

    buffer.set_size(nz, ny, nx);
    chasing_method_x->chasing(f, buffer); // Solve for x

    f.set_size(nz, nx, ny);
    buffer.transpose(f, {0, 2, 1}); // Transpose back

    buffer.set_size(nz, nx, ny);
    poisson_fft_y->transform_transpose(f, buffer);

    f.set_size(nx, ny, nz);
    buffer.transpose(f, {1, 2, 0});

    buffer.set_size(nx, ny, nz);
    poisson_fft_z->transform_transpose(f, buffer);

    std::swap(f, buffer);
    if (env_cfg.showCurrentStep)
        std::cout << "[Poisson] solve: done (domain " << domain_name << ")" << std::endl;
}

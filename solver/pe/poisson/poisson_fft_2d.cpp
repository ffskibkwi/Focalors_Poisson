#include "poisson_fft_2d.hpp"

void PoissonFFT2D<PDEBoundaryType::Dirichlet, PDEBoundaryType::Dirichlet>::init(int in_nx, int in_ny)
{
    nx = in_nx;
    ny = in_ny;

    sqr_2_N1 = std::sqrt(2.0 / (1. + ny));

    ex_vec = new double*[nx];
    for (int i = 0; i < nx; i++)
    {
        ex_vec[i] = new double[2 * ny + 2];
        for (int j = 0; j < 2 * ny + 2; j++)
            ex_vec[i][j] = 0.;
    }

    fft_result = new fftw_complex*[nx];
    for (int i = 0; i < nx; i++)
        fft_result[i] = new fftw_complex[ny + 2];

    plan = fftw_plan_dft_r2c_1d(2 * ny + 2, ex_vec[0], fft_result[0], FFTW_MEASURE);
}

PoissonFFT2D<PDEBoundaryType::Dirichlet, PDEBoundaryType::Dirichlet>::~PoissonFFT2D()
{
    fftw_destroy_plan(plan);
    for (int i = 0; i < nx; i++)
        delete[] ex_vec[i];
    delete[] ex_vec;
    for (int i = 0; i < nx; i++)
        delete[] fft_result[i];
    delete[] fft_result;
}

void PoissonFFT2D<PDEBoundaryType::Dirichlet, PDEBoundaryType::Dirichlet>::transform(const field2& f, field2& f_hat)
{
    OPENMP_PARALLEL_FOR()
    for (int i = 0; i < nx; i++)
    {
        for (int j = 0; j < ny; j++)
            ex_vec[i][j + 1] = f(i, j);
        fftw_execute_dft_r2c(plan, ex_vec[i], fft_result[i]);
        for (int j = 0; j < ny; j++)
            f_hat(i, j) = -sqr_2_N1 * fft_result[i][j + 1][1];
    }
}

void PoissonFFT2D<PDEBoundaryType::Dirichlet, PDEBoundaryType::Dirichlet>::transform_transpose(const field2& p_hat,
                                                                                                   field2&       p)
{
    OPENMP_PARALLEL_FOR()
    for (int i = 0; i < nx; i++)
    {
        for (int j = 0; j < ny; j++)
            ex_vec[i][j + 1] = p_hat(i, j);
        fftw_execute_dft_r2c(plan, ex_vec[i], fft_result[i]);
        for (int j = 0; j < ny; j++)
            p(i, j) = -sqr_2_N1 * fft_result[i][j + 1][1];
    }
}

void PoissonFFT2D<PDEBoundaryType::Neumann, PDEBoundaryType::Neumann>::init(int in_nx, int in_ny)
{
    nx = in_nx;
    ny = in_ny;

    sqr_2_n = std::sqrt(2.0 / ny);
    sqr_n   = std::sqrt(1.0 / ny);

    cos_ipi_n = new double[ny];
    sin_ipi_n = new double[ny];
    for (int j = 0; j < ny; j++)
    {
        cos_ipi_n[j] = std::cos(pi * (j + 1.) / ny / 2);
        sin_ipi_n[j] = std::sin(pi * (j + 1.) / ny / 2);
    }

    ex_vec = new double*[nx];
    for (int i = 0; i < nx; i++)
    {
        ex_vec[i] = new double[2 * ny];
        for (int j = 0; j < 2 * ny; j++)
            ex_vec[i][j] = 0.;
    }

    fft_result = new fftw_complex*[nx];
    for (int i = 0; i < nx; i++)
        fft_result[i] = new fftw_complex[ny + 1];

    plan = fftw_plan_dft_r2c_1d(2 * ny, ex_vec[0], fft_result[0], FFTW_MEASURE);
}

PoissonFFT2D<PDEBoundaryType::Neumann, PDEBoundaryType::Neumann>::~PoissonFFT2D()
{
    fftw_destroy_plan(plan);
    delete[] cos_ipi_n;
    delete[] sin_ipi_n;
    for (int i = 0; i < nx; i++)
        delete[] ex_vec[i];
    delete[] ex_vec;
    for (int i = 0; i < nx; i++)
        delete[] fft_result[i];
    delete[] fft_result;
}

void PoissonFFT2D<PDEBoundaryType::Neumann, PDEBoundaryType::Neumann>::transform(const field2& f, field2& f_hat)
{
    OPENMP_PARALLEL_FOR()
    for (int i = 0; i < nx; i++)
    {
        for (int j = 0; j < ny; j++)
            ex_vec[i][j + 1] = (j % 2 == 0 ? -1 : 1) * f(i, j);
        fftw_execute_dft_r2c(plan, ex_vec[i], fft_result[i]);
        for (int j = 0; j < (ny - 1); j++)
            f_hat(i, j) =
                sqr_2_n * (-1. * cos_ipi_n[j] * fft_result[i][j + 1][1] - sin_ipi_n[j] * fft_result[i][j + 1][0]);
        f_hat(i, ny - 1) = 0.;
        for (int j = 0; j < ny; j++)
            f_hat(i, ny - 1) += sqr_n * f(i, j);
    }
}

void PoissonFFT2D<PDEBoundaryType::Neumann, PDEBoundaryType::Neumann>::transform_transpose(const field2& p_hat,
                                                                                               field2&       p)
{
    OPENMP_PARALLEL_FOR()
    for (int i = 0; i < nx; i++)
    {
        for (int j = 0; j < ny; j++)
            ex_vec[i][j + 1] = cos_ipi_n[j] * p_hat(i, j);
        fftw_execute_dft_r2c(plan, ex_vec[i], fft_result[i]);
        for (int j = 0; j < ny; j++)
            p(i, j) = (j % 2 == 0 ? 1 : -1) * sqr_2_n * fft_result[i][j + 1][1];
        for (int j = 0; j < ny; j++)
            ex_vec[i][j + 1] = sin_ipi_n[j] * p_hat(i, j);
        fftw_execute_dft_r2c(plan, ex_vec[i], fft_result[i]);
        for (int j = 0; j < ny; j++)
            p(i, j) += (j % 2 == 0 ? 1 : -1) * sqr_2_n * fft_result[i][j + 1][0];
        for (int j = 0; j < ny; j++)
            p(i, j) += (sqr_2_n + sqr_n) * p_hat(i, ny - 1);
    }
}

void PoissonFFT2D<PDEBoundaryType::Neumann, PDEBoundaryType::Dirichlet>::init(int in_nx, int in_ny)
{
    nx = in_nx;
    ny = in_ny;

    sqr_4_2np1 = std::sqrt(4. / (2. * ny + 1.));

    cos_2nipi_2np1 = new double[ny];
    sin_2nipi_2np1 = new double[ny];
    for (int j = 0; j < ny; j++)
    {
        cos_2nipi_2np1[j] = std::cos(2. * ny * pi * (j + 1.) / (2. * ny + 1.));
        sin_2nipi_2np1[j] = std::sin(2. * ny * pi * (j + 1.) / (2. * ny + 1.));
    }

    ex_vec = new double*[nx];
    for (int i = 0; i < nx; i++)
    {
        ex_vec[i] = new double[2 * ny + 1];
        for (int j = 0; j < 2 * ny + 1; j++)
            ex_vec[i][j] = 0.;
    }

    fft_result = new fftw_complex*[nx];
    for (int i = 0; i < nx; i++)
        fft_result[i] = new fftw_complex[ny + 1];

    plan = fftw_plan_dft_r2c_1d(2 * ny + 1, ex_vec[0], fft_result[0], FFTW_MEASURE);
}

PoissonFFT2D<PDEBoundaryType::Neumann, PDEBoundaryType::Dirichlet>::~PoissonFFT2D()
{
    fftw_destroy_plan(plan);
    delete[] cos_2nipi_2np1;
    delete[] sin_2nipi_2np1;
    for (int i = 0; i < nx; i++)
        delete[] ex_vec[i];
    delete[] ex_vec;
    for (int i = 0; i < nx; i++)
        delete[] fft_result[i];
    delete[] fft_result;
}

void PoissonFFT2D<PDEBoundaryType::Neumann, PDEBoundaryType::Dirichlet>::transform(const field2& f, field2& f_hat)
{
    OPENMP_PARALLEL_FOR()
    for (int i = 0; i < nx; i++)
    {
        for (int j = 0; j < ny; j++)
            ex_vec[i][j + 1] = (j % 2 == 0 ? -1 : 1) * sqr_4_2np1 * f(i, j);
        fftw_execute_dft_r2c(plan, ex_vec[i], fft_result[i]);
        for (int j = 0; j < ny; j++)
            f_hat(i, j) = sin_2nipi_2np1[j] * fft_result[i][j + 1][0] - cos_2nipi_2np1[j] * fft_result[i][j + 1][1];
    }
}

void PoissonFFT2D<PDEBoundaryType::Neumann, PDEBoundaryType::Dirichlet>::transform_transpose(const field2& p_hat,
                                                                                                 field2&       p)
{
    OPENMP_PARALLEL_FOR()
    for (int i = 0; i < nx; i++)
    {
        for (int j = 0; j < ny; j++)
            ex_vec[i][j + 1] = sqr_4_2np1 * sin_2nipi_2np1[j] * p_hat(i, j);
        fftw_execute_dft_r2c(plan, ex_vec[i], fft_result[i]);
        for (int j = 0; j < ny; j++)
            p(i, j) = fft_result[i][j + 1][0];
        for (int j = 0; j < ny; j++)
            ex_vec[i][j + 1] = sqr_4_2np1 * cos_2nipi_2np1[j] * p_hat(i, j);
        fftw_execute_dft_r2c(plan, ex_vec[i], fft_result[i]);
        for (int j = 0; j < ny; j++)
            p(i, j) -= fft_result[i][j + 1][1];
        for (int j = 0; j < ny; j++)
            p(i, j) = (j % 2 == 0 ? -1 : 1) * p(i, j);
    }
}

void PoissonFFT2D<PDEBoundaryType::Dirichlet, PDEBoundaryType::Neumann>::init(int in_nx, int in_ny)
{
    nx = in_nx;
    ny = in_ny;

    sqr_4_2np1 = std::sqrt(4. / (2. * ny + 1.));

    ex_vec = new double*[nx];
    for (int i = 0; i < nx; i++)
    {
        ex_vec[i] = new double[2 * ny + 1];
        for (int j = 0; j < 2 * ny + 1; j++)
            ex_vec[i][j] = 0.;
    }

    fft_result = new fftw_complex*[nx];
    for (int i = 0; i < nx; i++)
        fft_result[i] = new fftw_complex[ny + 1];

    plan = fftw_plan_dft_r2c_1d(2 * ny + 1, ex_vec[0], fft_result[0], FFTW_MEASURE);
}

PoissonFFT2D<PDEBoundaryType::Dirichlet, PDEBoundaryType::Neumann>::~PoissonFFT2D()
{
    fftw_destroy_plan(plan);
    for (int i = 0; i < nx; i++)
        delete[] ex_vec[i];
    delete[] ex_vec;
    for (int i = 0; i < nx; i++)
        delete[] fft_result[i];
    delete[] fft_result;
}

void PoissonFFT2D<PDEBoundaryType::Dirichlet, PDEBoundaryType::Neumann>::transform(const field2& f, field2& f_hat)
{
    OPENMP_PARALLEL_FOR()
    for (int i = 0; i < nx; i++)
    {
        for (int j = 0; j < ny; j++)
            ex_vec[i][j + 1] = (j % 2 == 0 ? -1 : 1) * sqr_4_2np1 * f(i, j);
        fftw_execute_dft_r2c(plan, ex_vec[i], fft_result[i]);
        for (int j = 0; j < ny; j++)
            f_hat(i, j) = (j % 2 == 0 ? 1 : -1) * fft_result[i][j + 1][1];
    }
}

void PoissonFFT2D<PDEBoundaryType::Dirichlet, PDEBoundaryType::Neumann>::transform_transpose(const field2& p_hat,
                                                                                                 field2&       p)
{
    OPENMP_PARALLEL_FOR()
    for (int i = 0; i < nx; i++)
    {
        for (int j = 0; j < ny; j++)
            ex_vec[i][j + 1] = (j % 2 == 0 ? -1 : 1) * sqr_4_2np1 * p_hat(i, j);
        fftw_execute_dft_r2c(plan, ex_vec[i], fft_result[i]);
        for (int j = 0; j < ny; j++)
            p(i, j) = (j % 2 == 0 ? 1 : -1) * fft_result[i][j + 1][1];
    }
}

void PoissonFFT2D<PDEBoundaryType::Periodic, PDEBoundaryType::Periodic>::init(int in_nx, int in_ny)
{
    nx = in_nx;
    ny = in_ny;

    ex_vector_y_sin  = new double*[nx];
    ex_vector_y_cos  = new double*[nx];
    tc_y             = new double*[nx];
    ts_y             = new double*[nx];
    fft_result_y_sin = new fftw_complex*[nx];
    fft_result_y_cos = new fftw_complex*[nx];
    for (int i = 0; i < nx; i++)
    {
        ex_vector_y_sin[i]  = new double[ny];
        ex_vector_y_cos[i]  = new double[ny];
        tc_y[i]             = new double[ny];
        ts_y[i]             = new double[ny];
        fft_result_y_sin[i] = new fftw_complex[ny / 2 + 1];
        fft_result_y_cos[i] = new fftw_complex[ny / 2 + 1];
        for (int j = 0; j < ny; j++)
        {
            ex_vector_y_sin[i][j] = 0.0;
            ex_vector_y_cos[i][j] = 0.0;
            tc_y[i][j]            = 0.0;
            ts_y[i][j]            = 0.0;
        }
    }
    plan_y_sin = fftw_plan_dft_r2c_1d(ny, ex_vector_y_sin[0], fft_result_y_sin[0], FFTW_MEASURE);
    plan_y_cos = fftw_plan_dft_r2c_1d(ny, ex_vector_y_cos[0], fft_result_y_cos[0], FFTW_MEASURE);

    ex_vector_y_sin_even  = new double*[nx];
    ex_vector_y_cos_even  = new double*[nx];
    ex_vector_y_sin_odd   = new double*[nx];
    ex_vector_y_cos_odd   = new double*[nx];
    fc_y                  = new double*[nx];
    fs_y                  = new double*[nx];
    fft_result_y_sin_even = new fftw_complex*[nx];
    fft_result_y_cos_even = new fftw_complex*[nx];
    fft_result_y_sin_odd  = new fftw_complex*[nx];
    fft_result_y_cos_odd  = new fftw_complex*[nx];
    for (int i = 0; i < nx; i++)
    {
        ex_vector_y_sin_even[i]  = new double[ny];
        ex_vector_y_cos_even[i]  = new double[ny];
        ex_vector_y_sin_odd[i]   = new double[ny];
        ex_vector_y_cos_odd[i]   = new double[ny];
        fc_y[i]                  = new double[ny];
        fs_y[i]                  = new double[ny];
        fft_result_y_sin_even[i] = new fftw_complex[ny / 2 + 1];
        fft_result_y_cos_even[i] = new fftw_complex[ny / 2 + 1];
        fft_result_y_sin_odd[i]  = new fftw_complex[ny / 2 + 1];
        fft_result_y_cos_odd[i]  = new fftw_complex[ny / 2 + 1];
        for (int j = 0; j < ny; j++)
        {
            ex_vector_y_sin_even[i][j] = 0.0;
            ex_vector_y_cos_even[i][j] = 0.0;
            ex_vector_y_sin_odd[i][j]  = 0.0;
            ex_vector_y_cos_odd[i][j]  = 0.0;
            fc_y[i][j]                 = 0.0;
            fs_y[i][j]                 = 0.0;
        }
    }
    plan_y_sin_even = fftw_plan_dft_r2c_1d(ny, ex_vector_y_sin_even[0], fft_result_y_sin_even[0], FFTW_MEASURE);
    plan_y_cos_even = fftw_plan_dft_r2c_1d(ny, ex_vector_y_cos_even[0], fft_result_y_cos_even[0], FFTW_MEASURE);
    plan_y_sin_odd  = fftw_plan_dft_r2c_1d(ny, ex_vector_y_sin_odd[0], fft_result_y_sin_odd[0], FFTW_MEASURE);
    plan_y_cos_odd  = fftw_plan_dft_r2c_1d(ny, ex_vector_y_cos_odd[0], fft_result_y_cos_odd[0], FFTW_MEASURE);

    sqr_2_ny   = std::sqrt(2.0 / ny);
    cos_ipi_ny = new double[ny];
    sin_ipi_ny = new double[ny];
    for (int j = 0; j < ny; j++)
    {
        cos_ipi_ny[j] = std::cos(pi * (2.0 * j + 1.0) / 4.0);
        sin_ipi_ny[j] = std::sin(pi * (2.0 * j + 1.0) / 4.0);
    }
}

PoissonFFT2D<PDEBoundaryType::Periodic, PDEBoundaryType::Periodic>::~PoissonFFT2D()
{
    delete[] cos_ipi_ny;
    delete[] sin_ipi_ny;

    fftw_destroy_plan(plan_y_cos_odd);
    fftw_destroy_plan(plan_y_sin_odd);
    fftw_destroy_plan(plan_y_cos_even);
    fftw_destroy_plan(plan_y_sin_even);
    fftw_destroy_plan(plan_y_cos);
    fftw_destroy_plan(plan_y_sin);

    for (int i = 0; i < nx; i++)
    {
        delete[] ex_vector_y_sin_even[i];
        delete[] ex_vector_y_cos_even[i];
        delete[] ex_vector_y_sin_odd[i];
        delete[] ex_vector_y_cos_odd[i];
        delete[] fc_y[i];
        delete[] fs_y[i];
        delete[] fft_result_y_sin_even[i];
        delete[] fft_result_y_cos_even[i];
        delete[] fft_result_y_sin_odd[i];
        delete[] fft_result_y_cos_odd[i];
    }
    delete[] ex_vector_y_sin_even;
    delete[] ex_vector_y_cos_even;
    delete[] ex_vector_y_sin_odd;
    delete[] ex_vector_y_cos_odd;
    delete[] fc_y;
    delete[] fs_y;
    delete[] fft_result_y_sin_even;
    delete[] fft_result_y_cos_even;
    delete[] fft_result_y_sin_odd;
    delete[] fft_result_y_cos_odd;

    for (int i = 0; i < nx; i++)
    {
        delete[] ex_vector_y_sin[i];
        delete[] ex_vector_y_cos[i];
        delete[] tc_y[i];
        delete[] ts_y[i];
        delete[] fft_result_y_sin[i];
        delete[] fft_result_y_cos[i];
    }
    delete[] ex_vector_y_sin;
    delete[] ex_vector_y_cos;
    delete[] tc_y;
    delete[] ts_y;
    delete[] fft_result_y_sin;
    delete[] fft_result_y_cos;
}

void PoissonFFT2D<PDEBoundaryType::Periodic, PDEBoundaryType::Periodic>::transform(const field2& f, field2& f_hat)
\
{
    OPENMP_PARALLEL_FOR()
    for (int i = 0; i < nx; i++)
    {
        // For QT, only need one ex_vector_y
        for (int j = 0; j < ny - 1; j++)
        {
            ex_vector_y_sin[i][j + 1] = f(i, j);
        }
        fftw_execute_dft_r2c(plan_y_sin, ex_vector_y_sin[i], fft_result_y_sin[i]);
        for (int j = 0; j < ny / 2; j++)
        {
            if (j < ny / 2 - 1)
            {
                ts_y[i][2 * j + 2] = -1. * fft_result_y_sin[i][ny / 2 - j - 1][1];
            }
            ts_y[i][2 * j + 1] = -1. * fft_result_y_sin[i][ny / 2 - j - 1][1];
        }
        int curr_sign = 1;
        for (int j = 0; j < ny - 1; j++)
        {
            ex_vector_y_cos[i][j + 1] = curr_sign * f(i, j);
            curr_sign                 = -1 * curr_sign;
        }
        fftw_execute_dft_r2c(plan_y_cos, ex_vector_y_cos[i], fft_result_y_cos[i]);
        for (int j = 0; j < ny / 2; j++)
        {
            tc_y[i][2 * j]     = f(i, ny - 1) - fft_result_y_cos[i][j][0];
            tc_y[i][2 * j + 1] = f(i, ny - 1) - fft_result_y_cos[i][j + 1][0];
        }
        for (int j = 0; j < ny; j++)
        {
            f_hat(i, j) = sqr_2_ny * (cos_ipi_ny[j] * ts_y[i][j] + sin_ipi_ny[j] * tc_y[i][j]);
        }
    }
}

void PoissonFFT2D<PDEBoundaryType::Periodic, PDEBoundaryType::Periodic>::transform_transpose(const field2& p_hat,
                                                                                                 field2&       p)
\
{
    OPENMP_PARALLEL_FOR()
    for (int i = 0; i < nx; i++)
    {
        for (int j = 0; j < ny / 2; j++)
        {
            ex_vector_y_cos_even[i][j]             = cos_ipi_ny[2 * j] * p_hat(i, 2 * j);
            ex_vector_y_cos_odd[i][j + 1]          = cos_ipi_ny[2 * j + 1] * p_hat(i, 2 * j + 1);
            ex_vector_y_sin_even[i][ny / 2 - j]    = sin_ipi_ny[2 * j] * p_hat(i, 2 * j);
            ex_vector_y_sin_odd[i][ny / 2 - j - 1] = sin_ipi_ny[2 * j + 1] * p_hat(i, 2 * j + 1);
        }
        fftw_execute_dft_r2c(plan_y_cos_even, ex_vector_y_cos_even[i], fft_result_y_cos_even[i]);
        fftw_execute_dft_r2c(plan_y_cos_odd, ex_vector_y_cos_odd[i], fft_result_y_cos_odd[i]);
        fftw_execute_dft_r2c(plan_y_sin_even, ex_vector_y_sin_even[i], fft_result_y_sin_even[i]);
        fftw_execute_dft_r2c(plan_y_sin_odd, ex_vector_y_sin_odd[i], fft_result_y_sin_odd[i]);
        int curr_sign = -1;
        for (int j = 0; j < ny / 2; j++)
        {
            fc_y[i][j] = curr_sign * (fft_result_y_cos_even[i][j + 1][1] + fft_result_y_cos_odd[i][j + 1][1]);
            fs_y[i][j] = fft_result_y_sin_even[i][j + 1][0] + fft_result_y_sin_odd[i][j + 1][0];
            curr_sign  = -1 * curr_sign;
        }
        for (int j = 0; j < ny / 2 - 1; j++)
        {
            fc_y[i][(ny / 2) + j] = -1. * fc_y[i][ny / 2 - 2 - j];
            fs_y[i][(ny / 2) + j] = fs_y[i][ny / 2 - 2 - j];
        }
        fs_y[i][ny - 1] = fft_result_y_sin_even[i][0][0] + fft_result_y_sin_odd[i][0][0];
        for (int j = 0; j < ny; j++)
            p(i, j) = sqr_2_ny * (fc_y[i][j] + fs_y[i][j]);
    }
}
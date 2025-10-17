#include "poisson_fft_3d.hpp"

void PoissonFFT3D<PDEBoundaryType::Neumann, PDEBoundaryType::Neumann>::init(int in_nx, int in_ny, int in_nz)
{
    nx = in_nx;
    ny = in_ny;
    nz = in_nz;

    ex_vector_z  = new double*[nx];
    fft_result_z = new fftw_complex*[nx];
    for (int i = 0; i < nx; i++)
    {
        ex_vector_z[i]  = new double[nz * 2];
        fft_result_z[i] = new fftw_complex[nz * 2];
        for (int k = 0; k < nz * 2; ++k)
        {
            ex_vector_z[i][k] = 0.0;
        }
    }
    plan_z = fftw_plan_dft_r2c_1d(2 * nz, ex_vector_z[0], fft_result_z[0], FFTW_MEASURE);

    sqr_2_nz   = std::sqrt(2.0 / nz);
    sqr_nz     = std::sqrt(1.0 / nz);
    cos_ipi_nz = new double[nz];
    sin_ipi_nz = new double[nz];
    for (int k = 0; k < nz; k++)
    {
        cos_ipi_nz[k] = std::cos(pi * (k + 1) / nz / 2);
        sin_ipi_nz[k] = std::sin(pi * (k + 1) / nz / 2);
    }
}

PoissonFFT3D<PDEBoundaryType::Neumann, PDEBoundaryType::Neumann>::~PoissonFFT3D()
{
    delete[] cos_ipi_nz;
    delete[] sin_ipi_nz;

    for (int i = 0; i < nx; i++)
    {
        delete[] fft_result_z[i];
        delete[] ex_vector_z[i];
    }
    delete[] fft_result_z;
    delete[] ex_vector_z;
    fftw_destroy_plan(plan_z);
}

void PoissonFFT3D<PDEBoundaryType::Neumann, PDEBoundaryType::Neumann>::transform(const field3& f, field3& f_hat)
\
{
    OPENMP_PARALLEL_FOR()
    for (int i = 0; i < nx; i++)
    {
        for (int j = 0; j < ny; j++)
        {
            int curr_sign = -1;
            for (int k = 0; k < nz; k++)
            {
                ex_vector_z[i][k + 1] = curr_sign * f(i, j, k);
                // [k + 1]: The order of first num is 0, but in formula the calculation area is [1 , n]
                curr_sign = -1 * curr_sign;
            }
            fftw_execute_dft_r2c(plan_z, ex_vector_z[i], fft_result_z[i]);
            for (int k = 0; k < nz; k++)
            {
                f_hat(i, j, k) = sqr_2_nz * (-1. * cos_ipi_nz[k] * fft_result_z[i][k + 1][1] -
                                             sin_ipi_nz[k] * fft_result_z[i][k + 1][0]);
            }
            for (int k = 0; k < nz; k++)
                f_hat(i, j, nz - 1) += (sqr_nz + sqr_2_nz) * f(i, j, k);
        }
    }
}

void PoissonFFT3D<PDEBoundaryType::Neumann, PDEBoundaryType::Neumann>::transform_transpose(const field3& p_hat,
                                                                                               field3&       p)
\
{
    OPENMP_PARALLEL_FOR()
    for (int i = 0; i < nx; i++)
    {
        for (int j = 0; j < ny; j++)
        {
            for (int k = 0; k < nz; k++)
                ex_vector_z[i][k + 1] = cos_ipi_nz[k] * p_hat(i, j, k);
            fftw_execute_dft_r2c(plan_z, ex_vector_z[i], fft_result_z[i]);
            int curr_sign = 1;
            for (int k = 0; k < nz; k++)
            {
                p(i, j, k) = curr_sign * sqr_2_nz * fft_result_z[i][k + 1][1];
                curr_sign  = -1 * curr_sign;
            }
            for (int k = 0; k < nz; k++)
            {
                ex_vector_z[i][k + 1] = sin_ipi_nz[k] * p_hat(i, j, k);
            }
            fftw_execute_dft_r2c(plan_z, ex_vector_z[i], fft_result_z[i]);
            curr_sign = 1;
            for (int k = 0; k < nz; k++)
            {
                p(i, j, k) += curr_sign * sqr_2_nz * fft_result_z[i][k + 1][0];
                curr_sign = -1 * curr_sign;
            }
            for (int k = 0; k < nz; k++)
                p(i, j, k) += (sqr_2_nz + sqr_nz) * p_hat(i, j, nz - 1);
        }
    }
}

void PoissonFFT3D<PDEBoundaryType::Periodic, PDEBoundaryType::Periodic>::init(int in_nx, int in_ny, int in_nz)
{
    nx = in_nx;
    ny = in_ny;
    nz = in_nz;

    ex_vector_z_sin       = new double*[nx];
    ex_vector_z_cos       = new double*[nx];
    ex_vector_z_sin_even  = new double*[nx];
    ex_vector_z_cos_even  = new double*[nx];
    ex_vector_z_sin_odd   = new double*[nx];
    ex_vector_z_cos_odd   = new double*[nx];
    tc_z                  = new double*[nx];
    ts_z                  = new double*[nx];
    fc_z                  = new double*[nx];
    fs_z                  = new double*[nx];
    fft_result_z_sin      = new fftw_complex*[nx];
    fft_result_z_cos      = new fftw_complex*[nx];
    fft_result_z_sin_even = new fftw_complex*[nx];
    fft_result_z_cos_even = new fftw_complex*[nx];
    fft_result_z_sin_odd  = new fftw_complex*[nx];
    fft_result_z_cos_odd  = new fftw_complex*[nx];
    for (int i = 0; i < nx; i++)
    {
        ex_vector_z_sin[i]       = new double[nz];
        ex_vector_z_cos[i]       = new double[nz];
        ex_vector_z_sin_even[i]  = new double[nz];
        ex_vector_z_cos_even[i]  = new double[nz];
        ex_vector_z_sin_odd[i]   = new double[nz];
        ex_vector_z_cos_odd[i]   = new double[nz];
        tc_z[i]                  = new double[nz];
        ts_z[i]                  = new double[nz];
        fc_z[i]                  = new double[nz];
        fs_z[i]                  = new double[nz];
        fft_result_z_sin[i]      = new fftw_complex[nz / 2 + 1];
        fft_result_z_cos[i]      = new fftw_complex[nz / 2 + 1];
        fft_result_z_sin_even[i] = new fftw_complex[nz / 2 + 1];
        fft_result_z_cos_even[i] = new fftw_complex[nz / 2 + 1];
        fft_result_z_sin_odd[i]  = new fftw_complex[nz / 2 + 1];
        fft_result_z_cos_odd[i]  = new fftw_complex[nz / 2 + 1];
        for (int k = 0; k < nz; k++)
        {
            ex_vector_z_sin[i][k]      = 0.0;
            ex_vector_z_cos[i][k]      = 0.0;
            ex_vector_z_sin_even[i][k] = 0.0;
            ex_vector_z_cos_even[i][k] = 0.0;
            ex_vector_z_sin_odd[i][k]  = 0.0;
            ex_vector_z_cos_odd[i][k]  = 0.0;
            tc_z[i][k]                 = 0.0;
            ts_z[i][k]                 = 0.0;
            fc_z[i][k]                 = 0.0;
            fs_z[i][k]                 = 0.0;
        }
    }
    plan_z_sin      = fftw_plan_dft_r2c_1d(nz, ex_vector_z_sin[0], fft_result_z_sin[0], FFTW_MEASURE);
    plan_z_cos      = fftw_plan_dft_r2c_1d(nz, ex_vector_z_cos[0], fft_result_z_cos[0], FFTW_MEASURE);
    plan_z_sin_even = fftw_plan_dft_r2c_1d(nz, ex_vector_z_sin_even[0], fft_result_z_sin_even[0], FFTW_MEASURE);
    plan_z_cos_even = fftw_plan_dft_r2c_1d(nz, ex_vector_z_cos_even[0], fft_result_z_cos_even[0], FFTW_MEASURE);
    plan_z_sin_odd  = fftw_plan_dft_r2c_1d(nz, ex_vector_z_sin_odd[0], fft_result_z_sin_odd[0], FFTW_MEASURE);
    plan_z_cos_odd  = fftw_plan_dft_r2c_1d(nz, ex_vector_z_cos_odd[0], fft_result_z_cos_odd[0], FFTW_MEASURE);

    sqr_2_nz   = std::sqrt(2.0 / nz);
    sqr_nz     = std::sqrt(1.0 / nz);
    cos_ipi_nz = new double[nz];
    sin_ipi_nz = new double[nz];
    for (int k = 0; k < nz; k++)
    {
        cos_ipi_nz[k] = std::cos(pi * (2.0 * k + 1) / 4);
        sin_ipi_nz[k] = std::sin(pi * (2.0 * k + 1) / 4);
    }
}

PoissonFFT3D<PDEBoundaryType::Periodic, PDEBoundaryType::Periodic>::~PoissonFFT3D()
{
    delete[] sin_ipi_nz;
    delete[] cos_ipi_nz;

    for (int i = 0; i < nx; i++)
    {
        delete[] ex_vector_z_sin[i];
        delete[] ex_vector_z_cos[i];
        delete[] ex_vector_z_sin_even[i];
        delete[] ex_vector_z_cos_even[i];
        delete[] ex_vector_z_sin_odd[i];
        delete[] ex_vector_z_cos_odd[i];
        delete[] tc_z[i];
        delete[] ts_z[i];
        delete[] fc_z[i];
        delete[] fs_z[i];
        delete[] fft_result_z_sin[i];
        delete[] fft_result_z_cos[i];
        delete[] fft_result_z_sin_even[i];
        delete[] fft_result_z_cos_even[i];
        delete[] fft_result_z_sin_odd[i];
        delete[] fft_result_z_cos_odd[i];
    }
    delete[] ex_vector_z_sin;
    delete[] ex_vector_z_cos;
    delete[] ex_vector_z_sin_even;
    delete[] ex_vector_z_cos_even;
    delete[] ex_vector_z_sin_odd;
    delete[] ex_vector_z_cos_odd;
    delete[] tc_z;
    delete[] ts_z;
    delete[] fc_z;
    delete[] fs_z;
    delete[] fft_result_z_sin;
    delete[] fft_result_z_cos;
    delete[] fft_result_z_sin_even;
    delete[] fft_result_z_cos_even;
    delete[] fft_result_z_sin_odd;
    delete[] fft_result_z_cos_odd;
    fftw_destroy_plan(plan_z_sin);
    fftw_destroy_plan(plan_z_cos);
    fftw_destroy_plan(plan_z_sin_even);
    fftw_destroy_plan(plan_z_cos_even);
    fftw_destroy_plan(plan_z_sin_odd);
    fftw_destroy_plan(plan_z_cos_odd);
}

void PoissonFFT3D<PDEBoundaryType::Periodic, PDEBoundaryType::Periodic>::transform(const field3& f, field3& f_hat)
\
{
    OPENMP_PARALLEL_FOR()
    for (int i = 0; i < nx; i++)
    {
        for (int j = 0; j < ny; j++)
        {
            for (int k = 0; k < nz - 1; k++)
            {
                ex_vector_z_sin[i][k + 1] = f(i, j, k);
            }
            fftw_execute_dft_r2c(plan_z_sin, ex_vector_z_sin[i], fft_result_z_sin[i]);
            for (int k = 0; k < nz / 2; k++)
            {
                if (k < nz / 2 - 1)
                {
                    ts_z[i][2 * k + 2] = -1.0 * fft_result_z_sin[i][nz / 2 - k - 1][1];
                }
                ts_z[i][2 * k + 1] = -1.0 * fft_result_z_sin[i][nz / 2 - k - 1][1];
            }

            int curr_sign = 1;
            for (int k = 0; k < nz - 1; k++)
            {
                ex_vector_z_cos[i][k + 1] = curr_sign * f(i, j, k);
                curr_sign                 = -1 * curr_sign;
            }
            fftw_execute_dft_r2c(plan_z_cos, ex_vector_z_cos[i], fft_result_z_cos[i]);
            for (int k = 0; k < nz / 2; k++)
            {
                tc_z[i][2 * k]     = f(i, j, nz - 1) - fft_result_z_cos[i][k][0];
                tc_z[i][2 * k + 1] = f(i, j, nz - 1) - fft_result_z_cos[i][k + 1][0];
            }

            for (int k = 0; k < nz; k++)
            {
                f_hat(i, j, k) = sqr_2_nz * (cos_ipi_nz[k] * ts_z[i][k] + sin_ipi_nz[k] * tc_z[i][k]);
            }
        }
    }
}

void PoissonFFT3D<PDEBoundaryType::Periodic, PDEBoundaryType::Periodic>::transform_transpose(const field3& p_hat,
                                                                                                 field3&       p)
\
{
    OPENMP_PARALLEL_FOR()
    for (int i = 0; i < nx; i++)
    {
        for (int j = 0; j < ny; j++)
        {
            for (int k = 0; k < nz / 2; k++)
            {
                ex_vector_z_cos_even[i][k]             = cos_ipi_nz[2 * k] * p_hat(i, j, 2 * k);
                ex_vector_z_cos_odd[i][k + 1]          = cos_ipi_nz[2 * k + 1] * p_hat(i, j, 2 * k + 1);
                ex_vector_z_sin_even[i][nz / 2 - k]    = sin_ipi_nz[2 * k] * p_hat(i, j, 2 * k);
                ex_vector_z_sin_odd[i][nz / 2 - k - 1] = sin_ipi_nz[2 * k + 1] * p_hat(i, j, 2 * k + 1);
            }
            fftw_execute_dft_r2c(plan_z_cos_even, ex_vector_z_cos_even[i], fft_result_z_cos_even[i]);
            fftw_execute_dft_r2c(plan_z_cos_odd, ex_vector_z_cos_odd[i], fft_result_z_cos_odd[i]);
            fftw_execute_dft_r2c(plan_z_sin_even, ex_vector_z_sin_even[i], fft_result_z_sin_even[i]);
            fftw_execute_dft_r2c(plan_z_sin_odd, ex_vector_z_sin_odd[i], fft_result_z_sin_odd[i]);
            int curr_sign = -1;
            for (int k = 0; k < nz / 2; k++)
            {
                fc_z[i][k] = curr_sign * (fft_result_z_cos_even[i][k + 1][1] + fft_result_z_cos_odd[i][k + 1][1]);
                fs_z[i][k] = fft_result_z_sin_even[i][k + 1][0] + fft_result_z_sin_odd[i][k + 1][0];
                curr_sign  = -1 * curr_sign;
            }
            for (int k = 0; k < nz / 2 - 1; k++)
            {
                fc_z[i][nz / 2 + k] = -1. * fc_z[i][nz / 2 - 2 - k];
                fs_z[i][nz / 2 + k] = fs_z[i][nz / 2 - 2 - k];
            }
            fs_z[i][nz - 1] = fft_result_z_sin_even[i][0][0] + fft_result_z_sin_odd[i][0][0];
            for (int k = 0; k < nz; k++)
            {
                p(i, j, k) = sqr_2_nz * (fc_z[i][k] + fs_z[i][k]);
            }
        }
    }
}
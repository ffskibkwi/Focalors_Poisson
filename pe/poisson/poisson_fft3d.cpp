#include "poisson_fft3d.h"

// === Dirichlet-Dirichlet ===
PoissonFFT3D_DD::~PoissonFFT3D_DD()
{
    if (plan)
        fftw_destroy_plan(plan);
    if (ex_vec)
    {
        for (int i = 0; i < nx; ++i)
            delete[] ex_vec[i];
        delete[] ex_vec;
    }
    if (fft_result)
    {
        for (int i = 0; i < nx; ++i)
            delete[] fft_result[i];
        delete[] fft_result;
    }
}

void PoissonFFT3D_DD::init(int in_nx, int in_ny, int in_nz)
{
    nx       = in_nx;
    ny       = in_ny;
    nz       = in_nz;
    sqr_2_N1 = std::sqrt(2.0 / (1. + nz));

    ex_vec = new double*[nx];
    for (int i = 0; i < nx; i++)
    {
        ex_vec[i] = new double[2 * nz + 2];
        for (int k = 0; k < 2 * nz + 2; k++)
            ex_vec[i][k] = 0.;
    }

    fft_result = new fftw_complex*[nx];
    for (int i = 0; i < nx; i++)
        fft_result[i] = new fftw_complex[nz + 2];

    plan = fftw_plan_dft_r2c_1d(2 * nz + 2, ex_vec[0], fft_result[0], FFTW_MEASURE);
}

void PoissonFFT3D_DD::transform(const field3& f, field3& f_hat)
{
    OPENMP_PARALLEL_FOR()
    for (int i = 0; i < nx; i++)
    {
        for (int j = 0; j < ny; j++)
        {
            for (int k = 0; k < nz; k++)
                ex_vec[i][k + 1] = f(i, j, k);
            fftw_execute_dft_r2c(plan, ex_vec[i], fft_result[i]);
            for (int k = 0; k < nz; k++)
                f_hat(i, j, k) = -sqr_2_N1 * fft_result[i][k + 1][1];
        }
    }
}

void PoissonFFT3D_DD::transform_transpose(const field3& p_hat, field3& p)
{
    OPENMP_PARALLEL_FOR()
    for (int i = 0; i < nx; i++)
    {
        for (int j = 0; j < ny; j++)
        {
            for (int k = 0; k < nz; k++)
                ex_vec[i][k + 1] = p_hat(i, j, k);
            fftw_execute_dft_r2c(plan, ex_vec[i], fft_result[i]);
            for (int k = 0; k < nz; k++)
                p(i, j, k) = -sqr_2_N1 * fft_result[i][k + 1][1];
        }
    }
}

// === Neumann-Neumann ===
PoissonFFT3D_NN::~PoissonFFT3D_NN()
{
    if (plan)
        fftw_destroy_plan(plan);
    delete[] cos_ipi_n;
    delete[] sin_ipi_n;
    if (ex_vec)
    {
        for (int i = 0; i < nx; ++i)
            delete[] ex_vec[i];
        delete[] ex_vec;
    }
    if (fft_result)
    {
        for (int i = 0; i < nx; ++i)
            delete[] fft_result[i];
        delete[] fft_result;
    }
}

void PoissonFFT3D_NN::init(int in_nx, int in_ny, int in_nz)
{
    nx = in_nx;
    ny = in_ny;
    nz = in_nz;

    sqr_2_n = std::sqrt(2.0 / nz);
    sqr_n   = std::sqrt(1.0 / nz);

    cos_ipi_n = new double[nz];
    sin_ipi_n = new double[nz];
    for (int k = 0; k < nz; k++)
    {
        cos_ipi_n[k] = std::cos(pi * (k + 1.) / nz / 2);
        sin_ipi_n[k] = std::sin(pi * (k + 1.) / nz / 2);
    }

    ex_vec = new double*[nx];
    for (int i = 0; i < nx; i++)
    {
        ex_vec[i] = new double[2 * nz];
        for (int k = 0; k < 2 * nz; k++)
            ex_vec[i][k] = 0.;
    }

    fft_result = new fftw_complex*[nx];
    for (int i = 0; i < nx; i++)
        fft_result[i] = new fftw_complex[nz + 1];

    plan = fftw_plan_dft_r2c_1d(2 * nz, ex_vec[0], fft_result[0], FFTW_MEASURE);
}

void PoissonFFT3D_NN::transform(const field3& f, field3& f_hat)
{
    OPENMP_PARALLEL_FOR()
    for (int i = 0; i < nx; i++)
    {
        for (int j = 0; j < ny; j++)
        {
            for (int k = 0; k < nz; k++)
                ex_vec[i][k + 1] = (k % 2 == 0 ? -1 : 1) * f(i, j, k);
            fftw_execute_dft_r2c(plan, ex_vec[i], fft_result[i]);
            for (int k = 0; k < (nz - 1); k++)
                f_hat(i, j, k) =
                    sqr_2_n * (-1. * cos_ipi_n[k] * fft_result[i][k + 1][1] - sin_ipi_n[k] * fft_result[i][k + 1][0]);
            f_hat(i, j, nz - 1) = 0.;
            for (int k = 0; k < nz; k++)
                f_hat(i, j, nz - 1) += sqr_n * f(i, j, k);
        }
    }
}

void PoissonFFT3D_NN::transform_transpose(const field3& p_hat, field3& p)
{
    OPENMP_PARALLEL_FOR()
    for (int i = 0; i < nx; i++)
    {
        for (int j = 0; j < ny; j++)
        {
            for (int k = 0; k < nz; k++)
                ex_vec[i][k + 1] = cos_ipi_n[k] * p_hat(i, j, k);
            fftw_execute_dft_r2c(plan, ex_vec[i], fft_result[i]);
            for (int k = 0; k < nz; k++)
                p(i, j, k) = (k % 2 == 0 ? 1 : -1) * sqr_2_n * fft_result[i][k + 1][1];
            for (int k = 0; k < nz; k++)
                ex_vec[i][k + 1] = sin_ipi_n[k] * p_hat(i, j, k);
            fftw_execute_dft_r2c(plan, ex_vec[i], fft_result[i]);
            for (int k = 0; k < nz; k++)
                p(i, j, k) += (k % 2 == 0 ? 1 : -1) * sqr_2_n * fft_result[i][k + 1][0];
            for (int k = 0; k < nz; k++)
                p(i, j, k) += (sqr_2_n + sqr_n) * p_hat(i, j, nz - 1);
        }
    }
}

// === Neumann-Dirichlet ===
PoissonFFT3D_ND::~PoissonFFT3D_ND()
{
    if (plan)
        fftw_destroy_plan(plan);
    delete[] cos_2nipi_2np1;
    delete[] sin_2nipi_2np1;
    if (ex_vec)
    {
        for (int i = 0; i < nx; ++i)
            delete[] ex_vec[i];
        delete[] ex_vec;
    }
    if (fft_result)
    {
        for (int i = 0; i < nx; ++i)
            delete[] fft_result[i];
        delete[] fft_result;
    }
}

void PoissonFFT3D_ND::init(int in_nx, int in_ny, int in_nz)
{
    nx = in_nx;
    ny = in_ny;
    nz = in_nz;

    sqr_4_2np1 = std::sqrt(4. / (2. * nz + 1.));

    cos_2nipi_2np1 = new double[nz];
    sin_2nipi_2np1 = new double[nz];
    for (int k = 0; k < nz; k++)
    {
        cos_2nipi_2np1[k] = std::cos(2. * nz * pi * (k + 1.) / (2. * nz + 1.));
        sin_2nipi_2np1[k] = std::sin(2. * nz * pi * (k + 1.) / (2. * nz + 1.));
    }

    ex_vec = new double*[nx];
    for (int i = 0; i < nx; i++)
    {
        ex_vec[i] = new double[2 * nz + 1];
        for (int k = 0; k < 2 * nz + 1; k++)
            ex_vec[i][k] = 0.;
    }

    fft_result = new fftw_complex*[nx];
    for (int i = 0; i < nx; i++)
        fft_result[i] = new fftw_complex[nz + 1];

    plan = fftw_plan_dft_r2c_1d(2 * nz + 1, ex_vec[0], fft_result[0], FFTW_MEASURE);
}

void PoissonFFT3D_ND::transform(const field3& f, field3& f_hat)
{
    OPENMP_PARALLEL_FOR()
    for (int i = 0; i < nx; i++)
    {
        for (int j = 0; j < ny; j++)
        {
            for (int k = 0; k < nz; k++)
                ex_vec[i][k + 1] = (k % 2 == 0 ? -1 : 1) * sqr_4_2np1 * f(i, j, k);
            fftw_execute_dft_r2c(plan, ex_vec[i], fft_result[i]);
            for (int k = 0; k < nz; k++)
                f_hat(i, j, k) =
                    sin_2nipi_2np1[k] * fft_result[i][k + 1][0] - cos_2nipi_2np1[k] * fft_result[i][k + 1][1];
        }
    }
}

void PoissonFFT3D_ND::transform_transpose(const field3& p_hat, field3& p)
{
    OPENMP_PARALLEL_FOR()
    for (int i = 0; i < nx; i++)
    {
        for (int j = 0; j < ny; j++)
        {
            for (int k = 0; k < nz; k++)
                ex_vec[i][k + 1] = sqr_4_2np1 * sin_2nipi_2np1[k] * p_hat(i, j, k);
            fftw_execute_dft_r2c(plan, ex_vec[i], fft_result[i]);
            for (int k = 0; k < nz; k++)
                p(i, j, k) = fft_result[i][k + 1][0];
            for (int k = 0; k < nz; k++)
                ex_vec[i][k + 1] = sqr_4_2np1 * cos_2nipi_2np1[k] * p_hat(i, j, k);
            fftw_execute_dft_r2c(plan, ex_vec[i], fft_result[i]);
            for (int k = 0; k < nz; k++)
                p(i, j, k) -= fft_result[i][k + 1][1];
            for (int k = 0; k < nz; k++)
                p(i, j, k) = (k % 2 == 0 ? -1 : 1) * p(i, j, k);
        }
    }
}

// === Dirichlet-Neumann ===
PoissonFFT3D_DN::~PoissonFFT3D_DN()
{
    if (plan)
        fftw_destroy_plan(plan);
    if (ex_vec)
    {
        for (int i = 0; i < nx; ++i)
            delete[] ex_vec[i];
        delete[] ex_vec;
    }
    if (fft_result)
    {
        for (int i = 0; i < nx; ++i)
            delete[] fft_result[i];
        delete[] fft_result;
    }
}

void PoissonFFT3D_DN::init(int in_nx, int in_ny, int in_nz)
{
    nx = in_nx;
    ny = in_ny;
    nz = in_nz;

    sqr_4_2np1 = std::sqrt(4. / (2. * nz + 1.));

    ex_vec = new double*[nx];
    for (int i = 0; i < nx; i++)
    {
        ex_vec[i] = new double[2 * nz + 1];
        for (int k = 0; k < 2 * nz + 1; k++)
            ex_vec[i][k] = 0.;
    }

    fft_result = new fftw_complex*[nx];
    for (int i = 0; i < nx; i++)
        fft_result[i] = new fftw_complex[nz + 1];

    plan = fftw_plan_dft_r2c_1d(2 * nz + 1, ex_vec[0], fft_result[0], FFTW_MEASURE);
}

void PoissonFFT3D_DN::transform(const field3& f, field3& f_hat)
{
    OPENMP_PARALLEL_FOR()
    for (int i = 0; i < nx; i++)
    {
        for (int j = 0; j < ny; j++)
        {
            for (int k = 0; k < nz; k++)
                ex_vec[i][k + 1] = (k % 2 == 0 ? -1 : 1) * sqr_4_2np1 * f(i, j, k);
            fftw_execute_dft_r2c(plan, ex_vec[i], fft_result[i]);
            for (int k = 0; k < nz; k++)
                f_hat(i, j, k) = (k % 2 == 0 ? 1 : -1) * fft_result[i][k + 1][1];
        }
    }
}

void PoissonFFT3D_DN::transform_transpose(const field3& p_hat, field3& p)
{
    OPENMP_PARALLEL_FOR()
    for (int i = 0; i < nx; i++)
    {
        for (int j = 0; j < ny; j++)
        {
            for (int k = 0; k < nz; k++)
                ex_vec[i][k + 1] = (k % 2 == 0 ? -1 : 1) * sqr_4_2np1 * p_hat(i, j, k);
            fftw_execute_dft_r2c(plan, ex_vec[i], fft_result[i]);
            for (int k = 0; k < nz; k++)
                p(i, j, k) = (k % 2 == 0 ? 1 : -1) * fft_result[i][k + 1][1];
        }
    }
}

// === Periodic-Periodic ===
PoissonFFT3D_PP::~PoissonFFT3D_PP()
{
    delete[] cos_ipi_n;
    delete[] sin_ipi_n;

    if (plan_cos_odd)
        fftw_destroy_plan(plan_cos_odd);
    if (plan_sin_odd)
        fftw_destroy_plan(plan_sin_odd);
    if (plan_cos_even)
        fftw_destroy_plan(plan_cos_even);
    if (plan_sin_even)
        fftw_destroy_plan(plan_sin_even);
    if (plan_cos)
        fftw_destroy_plan(plan_cos);
    if (plan_sin)
        fftw_destroy_plan(plan_sin);

    if (ex_vector_sin_even)
    {
        for (int i = 0; i < nx; i++)
        {
            delete[] ex_vector_sin_even[i];
            delete[] ex_vector_cos_even[i];
            delete[] ex_vector_sin_odd[i];
            delete[] ex_vector_cos_odd[i];
            delete[] fc[i];
            delete[] fs[i];
            delete[] fft_result_sin_even[i];
            delete[] fft_result_cos_even[i];
            delete[] fft_result_sin_odd[i];
            delete[] fft_result_cos_odd[i];
        }
        delete[] ex_vector_sin_even;
        delete[] ex_vector_cos_even;
        delete[] ex_vector_sin_odd;
        delete[] ex_vector_cos_odd;
        delete[] fc;
        delete[] fs;
        delete[] fft_result_sin_even;
        delete[] fft_result_cos_even;
        delete[] fft_result_sin_odd;
        delete[] fft_result_cos_odd;
    }

    if (ex_vector_sin)
    {
        for (int i = 0; i < nx; i++)
        {
            delete[] ex_vector_sin[i];
            delete[] ex_vector_cos[i];
            delete[] tc[i];
            delete[] ts[i];
            delete[] fft_result_sin[i];
            delete[] fft_result_cos[i];
        }
        delete[] ex_vector_sin;
        delete[] ex_vector_cos;
        delete[] tc;
        delete[] ts;
        delete[] fft_result_sin;
        delete[] fft_result_cos;
    }
}

void PoissonFFT3D_PP::init(int in_nx, int in_ny, int in_nz)
{
    nx = in_nx;
    ny = in_ny;
    nz = in_nz;

    ex_vector_sin  = new double*[nx];
    ex_vector_cos  = new double*[nx];
    tc             = new double*[nx];
    ts             = new double*[nx];
    fft_result_sin = new fftw_complex*[nx];
    fft_result_cos = new fftw_complex*[nx];
    for (int i = 0; i < nx; i++)
    {
        ex_vector_sin[i]  = new double[nz];
        ex_vector_cos[i]  = new double[nz];
        tc[i]             = new double[nz];
        ts[i]             = new double[nz];
        fft_result_sin[i] = new fftw_complex[nz / 2 + 1];
        fft_result_cos[i] = new fftw_complex[nz / 2 + 1];
        for (int k = 0; k < nz; k++)
        {
            ex_vector_sin[i][k] = 0.0;
            ex_vector_cos[i][k] = 0.0;
            tc[i][k]            = 0.0;
            ts[i][k]            = 0.0;
        }
    }
    plan_sin = fftw_plan_dft_r2c_1d(nz, ex_vector_sin[0], fft_result_sin[0], FFTW_MEASURE);
    plan_cos = fftw_plan_dft_r2c_1d(nz, ex_vector_cos[0], fft_result_cos[0], FFTW_MEASURE);

    ex_vector_sin_even  = new double*[nx];
    ex_vector_cos_even  = new double*[nx];
    ex_vector_sin_odd   = new double*[nx];
    ex_vector_cos_odd   = new double*[nx];
    fc                  = new double*[nx];
    fs                  = new double*[nx];
    fft_result_sin_even = new fftw_complex*[nx];
    fft_result_cos_even = new fftw_complex*[nx];
    fft_result_sin_odd  = new fftw_complex*[nx];
    fft_result_cos_odd  = new fftw_complex*[nx];
    for (int i = 0; i < nx; i++)
    {
        ex_vector_sin_even[i]  = new double[nz];
        ex_vector_cos_even[i]  = new double[nz];
        ex_vector_sin_odd[i]   = new double[nz];
        ex_vector_cos_odd[i]   = new double[nz];
        fc[i]                  = new double[nz];
        fs[i]                  = new double[nz];
        fft_result_sin_even[i] = new fftw_complex[nz / 2 + 1];
        fft_result_cos_even[i] = new fftw_complex[nz / 2 + 1];
        fft_result_sin_odd[i]  = new fftw_complex[nz / 2 + 1];
        fft_result_cos_odd[i]  = new fftw_complex[nz / 2 + 1];
        for (int k = 0; k < nz; k++)
        {
            ex_vector_sin_even[i][k] = 0.0;
            ex_vector_cos_even[i][k] = 0.0;
            ex_vector_sin_odd[i][k]  = 0.0;
            ex_vector_cos_odd[i][k]  = 0.0;
            fc[i][k]                 = 0.0;
            fs[i][k]                 = 0.0;
        }
    }
    plan_sin_even = fftw_plan_dft_r2c_1d(nz, ex_vector_sin_even[0], fft_result_sin_even[0], FFTW_MEASURE);
    plan_cos_even = fftw_plan_dft_r2c_1d(nz, ex_vector_cos_even[0], fft_result_cos_even[0], FFTW_MEASURE);
    plan_sin_odd  = fftw_plan_dft_r2c_1d(nz, ex_vector_sin_odd[0], fft_result_sin_odd[0], FFTW_MEASURE);
    plan_cos_odd  = fftw_plan_dft_r2c_1d(nz, ex_vector_cos_odd[0], fft_result_cos_odd[0], FFTW_MEASURE);

    sqr_2_n   = std::sqrt(2.0 / nz);
    cos_ipi_n = new double[nz];
    sin_ipi_n = new double[nz];
    for (int k = 0; k < nz; k++)
    {
        cos_ipi_n[k] = std::cos(pi * (2.0 * k + 1.0) / 4.0);
        sin_ipi_n[k] = std::sin(pi * (2.0 * k + 1.0) / 4.0);
    }
}

void PoissonFFT3D_PP::transform(const field3& f, field3& f_hat)
{
    OPENMP_PARALLEL_FOR()
    for (int i = 0; i < nx; i++)
    {
        for (int j = 0; j < ny; j++)
        {
            for (int k = 0; k < nz - 1; k++)
                ex_vector_sin[i][k + 1] = f(i, j, k);
            fftw_execute_dft_r2c(plan_sin, ex_vector_sin[i], fft_result_sin[i]);
            for (int k = 0; k < nz / 2; k++)
            {
                if (k < nz / 2 - 1)
                {
                    ts[i][2 * k + 2] = -1. * fft_result_sin[i][nz / 2 - k - 1][1];
                }
                ts[i][2 * k + 1] = -1. * fft_result_sin[i][nz / 2 - k - 1][1];
            }
            int curr_sign = 1;
            for (int k = 0; k < nz - 1; k++)
            {
                ex_vector_cos[i][k + 1] = curr_sign * f(i, j, k);
                curr_sign               = -1 * curr_sign;
            }
            fftw_execute_dft_r2c(plan_cos, ex_vector_cos[i], fft_result_cos[i]);
            for (int k = 0; k < nz / 2; k++)
            {
                tc[i][2 * k]     = f(i, j, nz - 1) - fft_result_cos[i][k][0];
                tc[i][2 * k + 1] = f(i, j, nz - 1) - fft_result_cos[i][k + 1][0];
            }
            for (int k = 0; k < nz; k++)
            {
                f_hat(i, j, k) = sqr_2_n * (cos_ipi_n[k] * ts[i][k] + sin_ipi_n[k] * tc[i][k]);
            }
        }
    }
}

void PoissonFFT3D_PP::transform_transpose(const field3& p_hat, field3& p)
{
    OPENMP_PARALLEL_FOR()
    for (int i = 0; i < nx; i++)
    {
        for (int j = 0; j < ny; j++)
        {
            for (int k = 0; k < nz / 2; k++)
            {
                ex_vector_cos_even[i][k]             = cos_ipi_n[2 * k] * p_hat(i, j, 2 * k);
                ex_vector_cos_odd[i][k + 1]          = cos_ipi_n[2 * k + 1] * p_hat(i, j, 2 * k + 1);
                ex_vector_sin_even[i][nz / 2 - k]    = sin_ipi_n[2 * k] * p_hat(i, j, 2 * k);
                ex_vector_sin_odd[i][nz / 2 - k - 1] = sin_ipi_n[2 * k + 1] * p_hat(i, j, 2 * k + 1);
            }
            fftw_execute_dft_r2c(plan_cos_even, ex_vector_cos_even[i], fft_result_cos_even[i]);
            fftw_execute_dft_r2c(plan_cos_odd, ex_vector_cos_odd[i], fft_result_cos_odd[i]);
            fftw_execute_dft_r2c(plan_sin_even, ex_vector_sin_even[i], fft_result_sin_even[i]);
            fftw_execute_dft_r2c(plan_sin_odd, ex_vector_sin_odd[i], fft_result_sin_odd[i]);
            int curr_sign = -1;
            for (int k = 0; k < nz / 2; k++)
            {
                fc[i][k]  = curr_sign * (fft_result_cos_even[i][k + 1][1] + fft_result_cos_odd[i][k + 1][1]);
                fs[i][k]  = fft_result_sin_even[i][k + 1][0] + fft_result_sin_odd[i][k + 1][0];
                curr_sign = -1 * curr_sign;
            }
            for (int k = 0; k < nz / 2 - 1; k++)
            {
                fc[i][(nz / 2) + k] = -1. * fc[i][nz / 2 - 2 - k];
                fs[i][(nz / 2) + k] = fs[i][nz / 2 - 2 - k];
            }
            fs[i][nz - 1] = fft_result_sin_even[i][0][0] + fft_result_sin_odd[i][0][0];
            for (int k = 0; k < nz; k++)
                p(i, j, k) = sqr_2_n * (fc[i][k] + fs[i][k]);
        }
    }
}
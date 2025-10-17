#pragma once

#include "base/pch.h"

#include "base/location_boundary.h"
#include "poisson_fft_interface.h"

#include "fftw3.h"

/**
 * @brief fft in y direction on xy field.
 *
 * @tparam BoundTypeYNegative fluid boundary type in y negative.
 * @tparam BoundTypeYPositive fluid boundary type in y positive.
 */
template<PDEBoundaryType BoundTypeYNegative, PDEBoundaryType BoundTypeYPositive>
class PoissonFFT2D : public PoissonFFT2DInterface
{
public:
    PoissonFFT2D() {}

    void init(int in_nx, int in_ny) override
    {
        nx = in_nx;
        ny = in_ny;
    }

    void transform(const field2& f, field2& f_hat) override {}
    void transform_transpose(const field2& p_hat, field2& p) override {}

private:
    int nx, ny;
};

template<>
class PoissonFFT2D<PDEBoundaryType::Dirichlet, PDEBoundaryType::Dirichlet> : public PoissonFFT2DInterface
{
public:
    PoissonFFT2D() {}
    ~PoissonFFT2D();

    void init(int in_nx, int in_ny) override;
    void transform(const field2& f, field2& f_hat) override;
    void transform_transpose(const field2& p_hat, field2& p) override;

private:
    int nx = 0, ny = 0;

    double         sqr_2_N1   = 0.0;
    fftw_complex** fft_result = nullptr;
    fftw_plan      plan;
    double**       ex_vec = nullptr;
};

template<>
class PoissonFFT2D<PDEBoundaryType::Neumann, PDEBoundaryType::Neumann> : public PoissonFFT2DInterface
{
public:
    PoissonFFT2D() {}
    ~PoissonFFT2D();

    void init(int in_nx, int in_ny) override;
    void transform(const field2& f, field2& f_hat) override;
    void transform_transpose(const field2& p_hat, field2& p) override;

private:
    int nx = 0, ny = 0;

    double         sqr_2_n = 0.0, sqr_n = 0.0;
    double *       cos_ipi_n = nullptr, *sin_ipi_n = nullptr;
    fftw_complex** fft_result = nullptr;
    fftw_plan      plan;
    double**       ex_vec = nullptr;
};

template<>
class PoissonFFT2D<PDEBoundaryType::Neumann, PDEBoundaryType::Dirichlet> : public PoissonFFT2DInterface
{
public:
    PoissonFFT2D() {}
    ~PoissonFFT2D();

    void init(int in_nx, int in_ny) override;
    void transform(const field2& f, field2& f_hat) override;
    void transform_transpose(const field2& p_hat, field2& p) override;

private:
    int nx = 0, ny = 0;

    double         sqr_4_2np1     = 0.0;
    double *       cos_2nipi_2np1 = nullptr, *sin_2nipi_2np1 = nullptr;
    fftw_complex** fft_result = nullptr;
    fftw_plan      plan;
    double**       ex_vec = nullptr;
};

template<>
class PoissonFFT2D<PDEBoundaryType::Dirichlet, PDEBoundaryType::Neumann> : public PoissonFFT2DInterface
{
public:
    PoissonFFT2D() {}
    ~PoissonFFT2D();

    void init(int in_nx, int in_ny) override;
    void transform(const field2& f, field2& f_hat) override;
    void transform_transpose(const field2& p_hat, field2& p) override;

private:
    int nx = 0, ny = 0;

    double         sqr_4_2np1 = 0.0;
    fftw_complex** fft_result = nullptr;
    fftw_plan      plan;
    double**       ex_vec = nullptr;
};

template<>
class PoissonFFT2D<PDEBoundaryType::Periodic, PDEBoundaryType::Periodic> : public PoissonFFT2DInterface
{
public:
    PoissonFFT2D() {}
    ~PoissonFFT2D();

    void init(int in_nx, int in_ny) override;
    void transform(const field2& f, field2& f_hat) override;
    void transform_transpose(const field2& p_hat, field2& p) override;

private:
    int nx = 0, ny = 0;

    double**       ex_vector_y_sin  = nullptr;
    double**       ex_vector_y_cos  = nullptr;
    double**       tc_y             = nullptr;
    double**       ts_y             = nullptr;
    fftw_complex** fft_result_y_sin = nullptr;
    fftw_complex** fft_result_y_cos = nullptr;
    fftw_plan      plan_y_sin;
    fftw_plan      plan_y_cos;

    double**       ex_vector_y_sin_even  = nullptr;
    double**       ex_vector_y_cos_even  = nullptr;
    double**       ex_vector_y_sin_odd   = nullptr;
    double**       ex_vector_y_cos_odd   = nullptr;
    double**       fc_y                  = nullptr;
    double**       fs_y                  = nullptr;
    fftw_complex** fft_result_y_sin_even = nullptr;
    fftw_complex** fft_result_y_cos_even = nullptr;
    fftw_complex** fft_result_y_sin_odd  = nullptr;
    fftw_complex** fft_result_y_cos_odd  = nullptr;
    fftw_plan      plan_y_sin_even;
    fftw_plan      plan_y_cos_even;
    fftw_plan      plan_y_sin_odd;
    fftw_plan      plan_y_cos_odd;

    double* cos_ipi_nx = nullptr;
    double* sin_ipi_nx = nullptr;
    double* cos_ipi_ny = nullptr;
    double* sin_ipi_ny = nullptr;
    double  sqr_2_nx   = 0.0;
    double  sqr_2_ny   = 0.0;
};

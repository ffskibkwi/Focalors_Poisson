#pragma once

#include "pch.h"

#include "core/boundary/boundary_type.h"
#include "poisson_fft_interface.h"

#include "fftw3.h"

/**
 * @brief fft in z direction on xyz field.
 *
 * @tparam BoundTypeZNegative fluid boundary type in z negative.
 * @tparam BoundTypeZPositive fluid boundary type in z positive.
 */
template<PDEBoundaryType BoundTypeZNegative, PDEBoundaryType BoundTypeZPositive>
class PoissonFFT3D : public PoissonFFT3DInterface
{
public:
    PoissonFFT3D() {}

    void init(int in_nx, int in_ny, int in_nz) override
    {
        nx = in_nx;
        ny = in_ny;
        nz = in_nz;
    }

    void transform(const field3& f, field3& f_hat) override {}
    void transform_transpose(const field3& p_hat, field3& p) override {}

private:
    int nx, ny, nz;
};

template<>
class PoissonFFT3D<PDEBoundaryType::Neumann, PDEBoundaryType::Neumann> : public PoissonFFT3DInterface
{
public:
    PoissonFFT3D() {}
    ~PoissonFFT3D();

    void init(int in_nx, int in_ny, int in_nz) override;
    void transform(const field3& f, field3& f_hat) override;
    void transform_transpose(const field3& p_hat, field3& p) override;

private:
    int nx, ny, nz;

    double**       ex_vector_z  = nullptr;
    fftw_complex** fft_result_z = nullptr;
    fftw_plan      plan_z;

    double* cos_ipi_nz = nullptr; // cos(i*pi/nz)
    double* sin_ipi_nz = nullptr; // sin(i*pi/nz)
    double  sqr_2_nz   = 0.0;     // std::sqrt(2/nz)
    double  sqr_nz     = 0.0;     // std::sqrt(1/nz)
};

template<>
class PoissonFFT3D<PDEBoundaryType::Periodic, PDEBoundaryType::Periodic> : public PoissonFFT3DInterface
{
public:
    PoissonFFT3D() {}
    ~PoissonFFT3D();

    void init(int in_nx, int in_ny, int in_nz) override;
    void transform(const field3& f, field3& f_hat) override;
    void transform_transpose(const field3& p_hat, field3& p) override;

private:
    int nx, ny, nz;

    double**       ex_vector_z_sin       = nullptr;
    double**       ex_vector_z_cos       = nullptr;
    double**       ex_vector_z_sin_even  = nullptr;
    double**       ex_vector_z_cos_even  = nullptr;
    double**       ex_vector_z_sin_odd   = nullptr;
    double**       ex_vector_z_cos_odd   = nullptr;
    double**       tc_z                  = nullptr;
    double**       ts_z                  = nullptr;
    double**       fc_z                  = nullptr;
    double**       fs_z                  = nullptr;
    fftw_complex** fft_result_z_sin      = nullptr;
    fftw_complex** fft_result_z_cos      = nullptr;
    fftw_complex** fft_result_z_sin_even = nullptr;
    fftw_complex** fft_result_z_cos_even = nullptr;
    fftw_complex** fft_result_z_sin_odd  = nullptr;
    fftw_complex** fft_result_z_cos_odd  = nullptr;
    fftw_plan      plan_z_sin;
    fftw_plan      plan_z_cos;
    fftw_plan      plan_z_sin_even;
    fftw_plan      plan_z_cos_even;
    fftw_plan      plan_z_sin_odd;
    fftw_plan      plan_z_cos_odd;

    double* cos_ipi_nz = nullptr; // cos(i*pi/nz)
    double* sin_ipi_nz = nullptr; // sin(i*pi/nz)
    double  sqr_2_nz   = 0.0;     // std::sqrt(2/nz)
    double  sqr_nz     = 0.0;     // std::sqrt(1/nz)
};
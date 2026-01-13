#pragma once

#include "base/pch.h"

#include "base/location_boundary.h"
#include "fftw3.h"

class PoissonFFT3D
{
public:
    virtual ~PoissonFFT3D() {}

    virtual void init(int in_nx, int in_ny, int in_nz)               = 0;
    virtual void transform(const field3& f, field3& f_hat)           = 0;
    virtual void transform_transpose(const field3& p_hat, field3& p) = 0;
};

class PoissonFFT3D_DD : public PoissonFFT3D
{
public:
    PoissonFFT3D_DD() {}
    ~PoissonFFT3D_DD();

    void init(int in_nx, int in_ny, int in_nz) override;
    void transform(const field3& f, field3& f_hat) override;
    void transform_transpose(const field3& p_hat, field3& p) override;

private:
    int nx = 0, ny = 0, nz = 0;

    double         sqr_2_N1   = 0.0;
    fftw_complex** fft_result = nullptr;
    fftw_plan      plan       = nullptr;
    double**       ex_vec     = nullptr;
};

class PoissonFFT3D_NN : public PoissonFFT3D
{
public:
    PoissonFFT3D_NN() {}
    ~PoissonFFT3D_NN();

    void init(int in_nx, int in_ny, int in_nz) override;
    void transform(const field3& f, field3& f_hat) override;
    void transform_transpose(const field3& p_hat, field3& p) override;

private:
    int nx = 0, ny = 0, nz = 0;

    double         sqr_2_n = 0.0, sqr_n = 0.0;
    double *       cos_ipi_n = nullptr, *sin_ipi_n = nullptr;
    fftw_complex** fft_result = nullptr;
    fftw_plan      plan       = nullptr;
    double**       ex_vec     = nullptr;
};

class PoissonFFT3D_ND : public PoissonFFT3D
{
public:
    PoissonFFT3D_ND() {}
    ~PoissonFFT3D_ND();

    void init(int in_nx, int in_ny, int in_nz) override;
    void transform(const field3& f, field3& f_hat) override;
    void transform_transpose(const field3& p_hat, field3& p) override;

private:
    int nx = 0, ny = 0, nz = 0;

    double         sqr_4_2np1     = 0.0;
    double *       cos_2nipi_2np1 = nullptr, *sin_2nipi_2np1 = nullptr;
    fftw_complex** fft_result = nullptr;
    fftw_plan      plan       = nullptr;
    double**       ex_vec     = nullptr;
};

class PoissonFFT3D_DN : public PoissonFFT3D
{
public:
    PoissonFFT3D_DN() {}
    ~PoissonFFT3D_DN();

    void init(int in_nx, int in_ny, int in_nz) override;
    void transform(const field3& f, field3& f_hat) override;
    void transform_transpose(const field3& p_hat, field3& p) override;

private:
    int nx = 0, ny = 0, nz = 0;

    double         sqr_4_2np1 = 0.0;
    fftw_complex** fft_result = nullptr;
    fftw_plan      plan       = nullptr;
    double**       ex_vec     = nullptr;
};

class PoissonFFT3D_PP : public PoissonFFT3D
{
public:
    PoissonFFT3D_PP() {}
    ~PoissonFFT3D_PP();

    void init(int in_nx, int in_ny, int in_nz) override;
    void transform(const field3& f, field3& f_hat) override;
    void transform_transpose(const field3& p_hat, field3& p) override;

private:
    int nx = 0, ny = 0, nz = 0;

    double**       ex_vector_sin  = nullptr;
    double**       ex_vector_cos  = nullptr;
    double**       tc             = nullptr;
    double**       ts             = nullptr;
    fftw_complex** fft_result_sin = nullptr;
    fftw_complex** fft_result_cos = nullptr;
    fftw_plan      plan_sin       = nullptr;
    fftw_plan      plan_cos       = nullptr;

    double**       ex_vector_sin_even  = nullptr;
    double**       ex_vector_cos_even  = nullptr;
    double**       ex_vector_sin_odd   = nullptr;
    double**       ex_vector_cos_odd   = nullptr;
    double**       fc                  = nullptr;
    double**       fs                  = nullptr;
    fftw_complex** fft_result_sin_even = nullptr;
    fftw_complex** fft_result_cos_even = nullptr;
    fftw_complex** fft_result_sin_odd  = nullptr;
    fftw_complex** fft_result_cos_odd  = nullptr;
    fftw_plan      plan_sin_even       = nullptr;
    fftw_plan      plan_cos_even       = nullptr;
    fftw_plan      plan_sin_odd        = nullptr;
    fftw_plan      plan_cos_odd        = nullptr;

    double* cos_ipi_n = nullptr;
    double* sin_ipi_n = nullptr;
    double  sqr_2_n   = 0.0;
};

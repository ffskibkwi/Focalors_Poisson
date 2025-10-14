#pragma once

#include "pch.h"

#include "core/base/location_boundary.h"
#include "fftw3.h"

class PoissonFFT2D
{
public:
    virtual ~PoissonFFT2D() {}

    virtual void init(int in_nx, int in_ny)                          = 0;
    virtual void transform(const field2& f, field2& f_hat)           = 0;
    virtual void transform_transpose(const field2& p_hat, field2& p) = 0;
};

class PoissonFFT2D_DD : public PoissonFFT2D
{
public:
    PoissonFFT2D_DD() {}
    ~PoissonFFT2D_DD();

    void init(int in_nx, int in_ny) override;
    void transform(const field2& f, field2& f_hat) override;
    void transform_transpose(const field2& p_hat, field2& p) override;

private:
    int nx = 0, ny = 0;

    double         sqr_2_N1   = 0.0;
    fftw_complex** fft_result = nullptr;
    fftw_plan      plan       = nullptr;
    double**       ex_vec     = nullptr;
};

class PoissonFFT2D_NN : public PoissonFFT2D
{
public:
    PoissonFFT2D_NN() {}
    ~PoissonFFT2D_NN();

    void init(int in_nx, int in_ny) override;
    void transform(const field2& f, field2& f_hat) override;
    void transform_transpose(const field2& p_hat, field2& p) override;

private:
    int nx = 0, ny = 0;

    double         sqr_2_n = 0.0, sqr_n = 0.0;
    double *       cos_ipi_n = nullptr, *sin_ipi_n = nullptr;
    fftw_complex** fft_result = nullptr;
    fftw_plan      plan       = nullptr;
    double**       ex_vec     = nullptr;
};

class PoissonFFT2D_ND : public PoissonFFT2D
{
public:
    PoissonFFT2D_ND() {}
    ~PoissonFFT2D_ND();

    void init(int in_nx, int in_ny) override;
    void transform(const field2& f, field2& f_hat) override;
    void transform_transpose(const field2& p_hat, field2& p) override;

private:
    int nx = 0, ny = 0;

    double         sqr_4_2np1     = 0.0;
    double *       cos_2nipi_2np1 = nullptr, *sin_2nipi_2np1 = nullptr;
    fftw_complex** fft_result     = nullptr;
    fftw_plan      plan           = nullptr;
    double**       ex_vec         = nullptr;
};

class PoissonFFT2D_DN : public PoissonFFT2D
{
public:
    PoissonFFT2D_DN() {}
    ~PoissonFFT2D_DN();

    void init(int in_nx, int in_ny) override;
    void transform(const field2& f, field2& f_hat) override;
    void transform_transpose(const field2& p_hat, field2& p) override;

private:
    int nx = 0, ny = 0;

    double         sqr_4_2np1 = 0.0;
    fftw_complex** fft_result = nullptr;
    fftw_plan      plan       = nullptr;
    double**       ex_vec     = nullptr;
};

class PoissonFFT2D_PP : public PoissonFFT2D
{
public:
    PoissonFFT2D_PP() {}
    ~PoissonFFT2D_PP();

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
    fftw_plan      plan_y_sin       = nullptr;
    fftw_plan      plan_y_cos       = nullptr;

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
    fftw_plan      plan_y_sin_even       = nullptr;
    fftw_plan      plan_y_cos_even       = nullptr;
    fftw_plan      plan_y_sin_odd        = nullptr;
    fftw_plan      plan_y_cos_odd        = nullptr;

    double* cos_ipi_ny = nullptr;
    double* sin_ipi_ny = nullptr;
    double  sqr_2_ny   = 0.0;
};



#pragma once

#include "base/pch.h"

class PoissonFFT2DInterface
{
public:
    virtual void init(int in_nx, int in_ny)                          = 0;
    virtual void transform(const field2& f, field2& f_hat)           = 0;
    virtual void transform_transpose(const field2& p_hat, field2& p) = 0;
};

class PoissonFFT3DInterface
{
public:
    virtual void init(int in_nx, int in_ny, int in_nz)               = 0;
    virtual void transform(const field3& f, field3& f_hat)           = 0;
    virtual void transform_transpose(const field3& p_hat, field3& p) = 0;
};
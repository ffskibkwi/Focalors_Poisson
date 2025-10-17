#pragma once

#include "base/pch.h"

#include "output_scalable.h"

class PoissonSolver2DInterface : public OutputScalable
{
public:
    virtual void solve(field2& f) = 0;
};

class PoissonSolver3DInterface : public OutputScalable
{
public:
    virtual void solve(field3& f) = 0;
};
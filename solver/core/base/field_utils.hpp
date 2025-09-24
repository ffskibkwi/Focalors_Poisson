#pragma once

#include "pch.h"

#include <array>

namespace FieldUtils
{
    void transpose(const field2& src, field2& dst);
    void transpose(const field3& src, field3& dst, const std::array<int, 3>& permutation);
} // namespace FieldUtils

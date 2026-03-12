#pragma once

#include <cstdint>
#include <iostream>

enum class DifferenceSchemeType : std::uint8_t
{
    Conv_Center2nd_Diff_Center2nd,
    Conv_Upwind1st_Diff_Center2nd,
    Conv_QUICK_Diff_Center2nd,
    Conv_TVD_VanLeer_Diff_Center2nd,
};

std::ostream& operator<<(std::ostream& os, DifferenceSchemeType type);
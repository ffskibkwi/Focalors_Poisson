#include "scheme_type.h"

std::ostream& operator<<(std::ostream& os, DifferenceSchemeType type)
{
    switch (type)
    {
        case DifferenceSchemeType::Conv_Center2nd_Diff_Center2nd:
            os << "Conv_Center2nd_Diff_Center2nd";
            break;
        case DifferenceSchemeType::Conv_Upwind1st_Diff_Center2nd:
            os << "Conv_Upwind1st_Diff_Center2nd";
            break;
        case DifferenceSchemeType::Conv_QUICK_Diff_Center2nd:
            os << "Conv_QUICK_Diff_Center2nd";
            break;
        default:
            os << "Unknown";
            break;
    }
    return os;
}
#pragma once

#include "base/pch.h"

namespace IO
{
    bool array_to_csv(double** value, int nx, int ny, const std::string& filename);

    bool field_to_csv(field2& field, const std::string& filename);

    bool csv_to_field(field2& field, const std::string& filename);
} // namespace IO
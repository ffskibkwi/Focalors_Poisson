#pragma once

#include "base/pch.h"

namespace IO
{
    bool field_to_csv(field3& field, const std::string& filename);

    bool csv_to_field(field3& field, const std::string& filename);
} // namespace IO
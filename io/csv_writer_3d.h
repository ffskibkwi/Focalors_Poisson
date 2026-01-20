#pragma once

#include "base/pch.h"

namespace IO
{
    bool field_to_csv(field3& field, const std::string& filename);

    bool csv_to_field(field3& field, const std::string& filename);

    bool
    field_and_buffer_to_csv(field3& field, field2& buffer, const std::string& filename, VariablePositionType pos_type);

    bool var_to_csv(const Variable3D& var, const std::string& filename);
} // namespace IO
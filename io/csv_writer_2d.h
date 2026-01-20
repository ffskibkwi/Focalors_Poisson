#pragma once

#include "base/pch.h"

namespace IO
{
    bool array_to_csv(double** value, int nx, int ny, const std::string& filename);

    bool field_to_csv(field2& field, const std::string& filename);

    bool csv_to_field(field2& field, const std::string& filename);

    bool field_and_buffer_to_csv(field2& field, const std::string& filename, VariablePositionType pos_type);

    bool var_to_csv(const Variable& var, const std::string& filename);

    bool var_to_csv_full(const Variable& var, const std::string& filename);

    bool matlab_read_var(const Variable& var, const std::string& filename);
} // namespace IO

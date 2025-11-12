// paras_reader.h 新版参数读取接口
#pragma once

#include <string>
#include <unordered_map>

namespace IO
{
    std::unordered_map<std::string, std::string> paras_to_map(int argc, char* argv[]);
    bool read_string(const std::unordered_map<std::string, std::string>& para_map, const std::string& token, std::string& value_ref);
    bool read_bool(const std::unordered_map<std::string, std::string>& para_map, const std::string& token, bool& value_ref);

    template<typename T>
    bool read_number(const std::unordered_map<std::string, std::string>& para_map, const std::string& token, T& value_ref);
}
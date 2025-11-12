// paras_reader.cpp 新版参数读取实现
#include "paras_reader.h"
#include <algorithm>
#include <cctype>
#include <iostream>

namespace IO
{
    std::unordered_map<std::string, std::string> paras_to_map(int argc, char* argv[])
    {
        std::unordered_map<std::string, std::string> para_map;
        for (int i = 0; i < argc; ++i)
        {
            std::string line(argv[i]);
            int         index = line.find('=');
            if (index != std::string::npos)
            {
                para_map[line.substr(0, index)] = line.substr(index + 1);
            }
        }
        return para_map;
    }

    bool read_string(const std::unordered_map<std::string, std::string>& para_map,
                     const std::string&                                  token,
                     std::string&                                        value_ref)
    {
        auto it = para_map.find(token);
        if (it != para_map.end())
        {
            value_ref = it->second;
            return true;
        }
        return false;
    }

    bool
    read_bool(const std::unordered_map<std::string, std::string>& para_map, const std::string& token, bool& value_ref)
    {
        std::string value_str;
        if (read_string(para_map, token, value_str))
        {
            std::transform(
                value_str.begin(), value_str.end(), value_str.begin(), [](unsigned char c) { return std::tolower(c); });
            if (value_str == "true" || value_str == "1")
            {
                value_ref = true;
                return true;
            }
            else if (value_str == "false" || value_str == "0")
            {
                value_ref = false;
                return true;
            }
            else
            {
                std::cerr << "Error: Invalid boolean value for token '" << token << "'." << std::endl;
                return false;
            }
        }
        return false;
    }

    template<typename T>
    bool
    read_number(const std::unordered_map<std::string, std::string>& para_map, const std::string& token, T& value_ref)
    {
        auto it = para_map.find(token);
        if (it != para_map.end())
        {
            const std::string& value_str = it->second;
            size_t             index     = value_str.find('/');
            try
            {
                if (index != std::string::npos)
                {
                    double numerator   = std::stod(value_str.substr(0, index));
                    double denominator = std::stod(value_str.substr(index + 1));
                    if (denominator == 0)
                    {
                        std::cerr << "Error: Division by zero." << std::endl;
                        return false;
                    }
                    value_ref = static_cast<T>(numerator / denominator);
                }
                else
                {
                    value_ref = static_cast<T>(std::stod(value_str));
                }
                return true;
            }
            catch (const std::exception& e)
            {
                std::cerr << "Error: Invalid argument for token '" << token << "'." << std::endl;
                return false;
            }
        }
        return false;
    }
} // namespace IO
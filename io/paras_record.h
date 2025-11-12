// paras_record.h 新版参数记录接口
#pragma once

#include <fstream>
#include <string>

namespace IO
{
    enum SimInfoOutputFormat
    {
        Matlab,
        Python,
        Default
    };

    class ParasRecord
    {
    public:
        ParasRecord() = default;
        ParasRecord(const std::string& root_dir, SimInfoOutputFormat output_format);
        ~ParasRecord();

        template<typename ValueType>
        ParasRecord& record(const std::string& key, ValueType value);

        ParasRecord& record(const std::string& key, const std::string& value);

    private:
        std::string   output_filename;
        std::ofstream outfile;
        std::string   delimiter = ": ";
        std::string   terminator;
    };
} // namespace IO
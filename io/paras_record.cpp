// paras_record.cpp 新版参数记录实现
#include "paras_record.h"
#include <sstream>

namespace IO
{
    ParasRecord::ParasRecord(const std::string& root_dir, SimInfoOutputFormat output_format)
    {
        std::stringstream ss;
        ss << root_dir;
        std::string extend;
        switch (output_format)
        {
            case SimInfoOutputFormat::Matlab:
                extend = ".m";
                break;
            case SimInfoOutputFormat::Python:
                extend = ".py";
                break;
            default:
                extend = ".txt";
                break;
        }
        ss << "/paras_record" << extend;
        output_filename = ss.str();
        outfile.open(output_filename, std::ios::app);

        switch (output_format)
        {
            case SimInfoOutputFormat::Matlab:
                delimiter  = " = ";
                terminator = ";";
                break;
            case SimInfoOutputFormat::Python:
                delimiter  = " = ";
                terminator = "";
                break;
            default:
                delimiter  = ": ";
                terminator = "";
                break;
        }
    }

    ParasRecord::~ParasRecord()
    {
        if (outfile.is_open())
            outfile.close();
    }

    template<typename ValueType>
    ParasRecord& ParasRecord::record(const std::string& key, ValueType value)
    {
        outfile << key << delimiter << value << terminator << std::endl;
        return *this;
    }

    ParasRecord& ParasRecord::record(const std::string& key, const std::string& value)
    {
        outfile << key << delimiter << "\"" << value << "\"" << terminator << std::endl;
        return *this;
    }
} // namespace IO
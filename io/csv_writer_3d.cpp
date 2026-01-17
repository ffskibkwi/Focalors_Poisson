#include "csv_writer_3d.h"
#include "common.h"

#include <fstream>

namespace IO
{
    bool field_to_csv(field3& field, const std::string& filename)
    {
        fs::path path(filename);
        fs::path dir = path.parent_path();
        IO::create_directory(dir);

        std::ofstream outfile(filename + ".csv");

        if (!outfile.is_open())
        {
            std::cerr << "Failed to open file: " << filename + ".csv" << std::endl;
            return false;
        }

        for (int i = 0; i < field.get_nx(); i++)
        {
            for (int j = 0; j < field.get_ny(); j++)
            {
                for (int k = 0; k < field.get_nz(); k++)
                {
                    if (k == (field.get_nz() - 1))
                    {
                        outfile << field(i, j, k) << std::endl;
                    }
                    else
                    {
                        outfile << field(i, j, k) << ",";
                    }
                }
            }
        }
        outfile.close();

        if (outfile.fail())
        {
            std::cout << "Writing to file failed: " << std::endl;
            return false;
        }

        return true;
    }

    bool csv_to_field(field3& field, const std::string& filename)
    {
        std::ifstream infile(filename + ".csv");

        if (!infile.is_open())
        {
            std::cerr << "Failed to open file: " << filename + ".csv" << std::endl;
            return false;
        }

        std::string line;
        int         i = 0, j = 0;

        while (std::getline(infile, line))
        {
            std::stringstream ss(line);
            std::string       value;
            int               k = 0;

            while (std::getline(ss, value, ','))
            {
                try
                {
                    double numeric_value = std::stod(value);
                    field(i, j, k)       = numeric_value;
                    k++;
                }
                catch (const std::invalid_argument&)
                {
                    std::cerr << "Invalid number at i = " << i << ", j = " << j << ", k = " << k << ": " << value
                              << std::endl;
                }
            }
            j++;
            if (j >= field.get_ny())
            {
                j = 0;
                i++;
            }
        }

        infile.close();
        return true;
    }
} // namespace IO
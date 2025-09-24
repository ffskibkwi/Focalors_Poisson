#include "csv_writer_2d.h"
#include "common.h"

#include <fstream>

namespace IO
{
    bool array_to_csv(double** value, int nx, int ny, const std::string& filename)
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

        for (int i = 0; i < nx; i++)
        {
            for (int j = 0; j < ny; j++)
            {
                if (j == (ny - 1))
                {
                    outfile << value[i][j] << std::endl;
                }
                else
                {
                    outfile << value[i][j] << ",";
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

    bool field_to_csv(field2& field, const std::string& filename)
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
                if (j == (field.get_ny() - 1))
                {
                    outfile << field(i, j) << std::endl;
                }
                else
                {
                    outfile << field(i, j) << ",";
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

    bool csv_to_field(field2& field, const std::string& filename)
    {
        std::ifstream infile(filename + ".csv");

        if (!infile.is_open())
        {
            std::cerr << "Failed to open file: " << filename + ".csv" << std::endl;
            return false;
        }

        std::string line;
        int         i = 0;

        while (std::getline(infile, line))
        {
            std::stringstream ss(line);
            std::string       value;
            int               j = 0;

            while (std::getline(ss, value, ','))
            {
                try
                {
                    double numeric_value = std::stod(value);
                    field(i, j)          = numeric_value;
                    j++;
                }
                catch (const std::invalid_argument&)
                {
                    std::cerr << "Invalid number at i " << i << ", j " << j << ": " << value << std::endl;
                }
            }
            i++;
        }

        infile.close();
        return true;
    }
} // namespace IO
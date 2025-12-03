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

    bool
    field_and_buffer_to_csv(field2& field, double* buffer, const std::string& filename, VariablePositionType pos_type)
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

        int nx = field.get_nx();
        int ny = field.get_ny();

        if (pos_type == VariablePositionType::XEdge)
        {
            for (int i = 0; i < nx + 1; i++)
            {
                for (int j = 0; j < ny; j++)
                {
                    double val = 0;
                    if (i == nx)
                    {
                        val = buffer[j];
                    }
                    else
                    {
                        val = field(i, j);
                    }

                    if (j == (ny - 1))
                    {
                        outfile << val << std::endl;
                    }
                    else
                    {
                        outfile << val << ",";
                    }
                }
            }
        }
        else if (pos_type == VariablePositionType::YEdge)
        {
            for (int i = 0; i < nx; i++)
            {
                for (int j = 0; j < ny + 1; j++)
                {
                    double val = 0;
                    if (j == ny)
                    {
                        val = buffer[i];
                    }
                    else
                    {
                        val = field(i, j);
                    }

                    if (j == ny)
                    {
                        outfile << val << std::endl;
                    }
                    else
                    {
                        outfile << val << ",";
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

    bool var_to_csv(const Variable& var, const std::string& filename)
    {
        auto& domains        = var.geometry->domains;
        auto& field_map      = var.field_map;
        auto& buffer_map     = var.buffer_map;
        auto& boundary_types = var.boundary_type_map;

        for (auto& domain : domains)
        {
            try
            {
                auto& field         = field_map.at(domain);
                auto& buffers       = buffer_map.at(domain);
                auto& boundary_type = boundary_types.at(domain);

                int nx = field->get_nx();
                int ny = field->get_ny();
                if (var.position_type == VariablePositionType::XEdge)
                {
                    if (boundary_type.at(LocationType::Right) == PDEBoundaryType::Adjacented)
                        field_to_csv(*field, filename + "_" + domain->name);
                    else
                        field_and_buffer_to_csv(*field,
                                                buffers.at(LocationType::Right),
                                                filename + "_" + domain->name,
                                                VariablePositionType::XEdge);
                }
                else if (var.position_type == VariablePositionType::YEdge)
                {
                    if (boundary_type.at(LocationType::Up) == PDEBoundaryType::Adjacented)
                        field_to_csv(*field, filename + "_" + domain->name);
                    else
                        field_and_buffer_to_csv(*field,
                                                buffers.at(LocationType::Up),
                                                filename + "_" + domain->name,
                                                VariablePositionType::YEdge);
                }
                else
                {
                    field_to_csv(*field, filename + "_" + domain->name);
                }
            }
            catch (const std::exception& e)
            {
                std::cerr << "[var_to_csv] Error: " << e.what() << std::endl;
                return false;
            }
        }
        return true;
    }

    bool var_to_csv_full(const Variable& var, const std::string& filename)
    {
        auto& domains    = var.geometry->domains;
        auto& field_map  = var.field_map;
        auto& buffer_map = var.buffer_map;

        for (auto& domain : domains)
        {
            try
            {
                auto& field   = field_map.at(domain);
                auto& buffers = buffer_map.at(domain);

                if (var.position_type == VariablePositionType::XEdge)
                {
                    field_and_buffer_to_csv(*field,
                                            buffers.at(LocationType::Right),
                                            filename + "_" + domain->name,
                                            VariablePositionType::XEdge);
                }
                else if (var.position_type == VariablePositionType::YEdge)
                {
                    field_and_buffer_to_csv(*field,
                                            buffers.at(LocationType::Up),
                                            filename + "_" + domain->name,
                                            VariablePositionType::YEdge);
                }
                else
                {
                    field_to_csv(*field, filename + "_" + domain->name);
                }
            }
            catch (const std::exception& e)
            {
                std::cerr << "[var_to_csv_full] Error: " << e.what() << std::endl;
                return false;
            }
        }
        return true;
    }
} // namespace IO
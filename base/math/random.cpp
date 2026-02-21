#include "random.h"

void add_random_number(double* field, int length, double min, double max, int seed)
{
    std::mt19937 gen;
    if (seed == -1)
    {
        gen = std::mt19937(std::random_device {}());
    }
    else
    {
        gen = std::mt19937(seed);
    }

    std::uniform_real_distribution<> dis(min, max);

    // no opemp for reproducible random numbers
    for (int i = 0; i < length; i++)
    {
        field[i] += dis(gen);
    }
}

void add_random_number(field2& field, double min, double max, int seed)
{
    std::mt19937 gen;
    if (seed == -1)
    {
        gen = std::mt19937(std::random_device {}());
    }
    else
    {
        gen = std::mt19937(seed);
    }

    std::uniform_real_distribution<> dis(min, max);

    // no opemp for reproducible random numbers
    for (int i = 1; i < field.get_nx() - 1; i++)
    {
        for (int j = 1; j < field.get_ny() - 1; j++)
        {
            field(i, j) += dis(gen);
        }
    }
}

void add_random_number(field3& field, double min, double max, int seed)
{
    std::mt19937 gen;
    if (seed == -1)
    {
        gen = std::mt19937(std::random_device {}());
    }
    else
    {
        gen = std::mt19937(seed);
    }

    std::uniform_real_distribution<> dis(min, max);

    // no opemp for reproducible random numbers
    for (int i = 1; i < field.get_nx() - 1; i++)
    {
        for (int j = 1; j < field.get_ny() - 1; j++)
        {
            for (int k = 1; k < field.get_nz() - 1; k++)
            {
                field(i, j, k) += dis(gen);
            }
        }
    }
}
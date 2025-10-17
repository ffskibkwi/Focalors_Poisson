#pragma once

#include <string>
#include "base/location_boundary.h"

class field2
{
public:
    double* value = nullptr;

    field2() {}
    field2(int in_nx, int in_ny);
    field2(int in_nx, int in_ny, const std::string& in_name);
    field2(const std::string& in_name);
    ~field2();
    field2(const field2& rhs) noexcept;
    field2& operator=(const field2& rhs) noexcept;
    field2(field2&& rhs) noexcept { swap(*this, rhs); }
    field2& operator=(field2&& rhs) noexcept
    {
        if (this != &rhs)
        {
            swap(*this, rhs);
        }
        return *this;
    }

    void init(int in_nx, int in_ny, const std::string& in_name);
    void init(int in_nx, int in_ny);
    void deinit();

    // Basic calculation
    field2 operator+(const field2& rhs);
    field2 operator-(const field2& rhs);
    field2 operator*(const double k);
    field2 operator*(const field2& rhs);
    void   add(const double t);
    void   add_affine_transform(const double a, const field2& x, const double b);
    double dot(const field2& A);
    double norm();
    double sum();
    double squared_sum();
    double mean_at_x_axis(int i);
    double mean_at_y_axis(int j);

    // Boundary operation
    void left_bond_add(const double k, const field2& left_field2);
    void right_bond_add(const double k, const field2& right_field2);
    void up_bond_add(const double k, const field2& up_field2);
    void down_bond_add(const double k, const field2& down_field2);
    void bond_add(LocationType location, const double k, const field2& neighbour_field2);

    double& operator()(int i, int j);
    double  operator()(int i, int j) const;

    int               get_nx() const { return nx; }
    int               get_ny() const { return ny; }
    int               get_size_n() const { return size_n; }
    bool              set_size(int in_nx, int in_ny);
    double*           get_ptr(int i, unsigned int j);
    double*           get_ptr(int i);
    const std::string get_name() const { return name; }

    void clear(double clear_value = 0.0);

    friend void swap(field2& lhs, field2& rhs);
    void transpose(field2& dst);

protected:
    unsigned int nx, ny;
    unsigned int size_n;
    std::string  name = "Default";
};
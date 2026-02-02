#pragma once

#include "base/location_boundary.h"

#include <array>
#include <string>

class field2;
class field3
{
public:
    double* value = nullptr;

    field3() {}
    field3(int in_nx, int in_ny, int in_nz);
    field3(int in_nx, int in_ny, int in_nz, const std::string& in_name);
    field3(const std::string& in_name);
    ~field3();
    field3& operator=(const field3& rhs);
    field3(field3&& rhs) noexcept { swap(*this, rhs); }
    field3& operator=(field3&& rhs) noexcept
    {
        if (this != &rhs)
        {
            swap(*this, rhs);
        }
        return *this;
    }

    void init(int in_nx, int in_ny, int in_nz, const std::string& in_name = "");
    void deinit();

    field3  operator+(const field3& rhs);
    field3  operator-(const field3& rhs);
    field3  operator*(const double a);
    field3& operator*=(const double a);
    void    add_affine_transform(const double a, const field3& x, const double b);
    void    transpose(field3& dst, const std::array<int, 3>& permutation);

    double norm();

    double sum();
    // Computes the sum of elements in a xy-plane at a given z-index.
    double sum_at_xy_plane(int k);
    // Computes the sum of elements in a xz-plane at a given y-index.
    double sum_at_xz_plane(int j);
    // Computes the sum of elements in a yz-plane at a given x-index.
    double sum_at_yz_plane(int i);

    double mean();
    // Computes the mean of elements in a xy-plane at a given z-index.
    double mean_at_xy_plane(int k);
    // Computes the mean of elements in a xz-plane at a given y-index.
    double mean_at_xz_plane(int j);
    // Computes the mean of elements in a yz-plane at a given x-index.
    double mean_at_yz_plane(int i);

    // Boundary operation
    void left_bond_add(const double a, const field2& bound);
    void right_bond_add(const double a, const field2& bound);
    void front_bond_add(const double a, const field2& bound);
    void back_bond_add(const double a, const field2& bound);
    void down_bond_add(const double a, const field2& bound);
    void up_bond_add(const double a, const field2& bound);
    void bond_add(LocationType location, const double a, const field2& bound);

    void left_bond_add(const double a, const field3& neighbour);
    void right_bond_add(const double a, const field3& neighbour);
    void front_bond_add(const double a, const field3& neighbour);
    void back_bond_add(const double a, const field3& neighbour);
    void down_bond_add(const double a, const field3& neighbour);
    void up_bond_add(const double a, const field3& neighbour);
    void bond_add(LocationType location, const double a, const field3& neighbour);

    double& operator()(int i, int j, int k);
    double  operator()(int i, int j, int k) const;

    int               get_nx() const { return nx; }
    int               get_ny() const { return ny; }
    int               get_nz() const { return nz; }
    int               get_size_n() const { return size_n; }
    bool              set_size(int in_nx, int in_ny, int in_nz);
    double*           get_ptr(int i, int j, int k) const;
    const std::string get_name() const { return name; }

    /**
     * @brief Copies a slice of xz-plane data to a destination array.
     *
     * @param dest      Destination array for the copied data.
     * @param start     Starting y-index of the slice.
     * @param length    Number of y-planes to copy.
     */
    void copy_xz_slice_to(double* dest, int start, int length);

    /**
     * @brief Copies data from a source array to a slice of xz-planes.
     *
     * @param src       Source array containing the data to copy.
     * @param start     Starting y-index of the slice.
     * @param length    Number of y-planes to copy.
     */
    void copy_xz_slice_from(double* src, int start, int length);

    /**
     * @brief Copies a block of xy-plane data to a destination array.
     *
     * @param dest      Destination array for the copied data.
     * @param start_x   Starting x-index of the block.
     * @param start_y   Starting y-index of the block.
     * @param length_x  Number of x-planes to copy.
     * @param length_y  Number of y-planes to copy.
     */
    void copy_xy_block_to(double* dest, int start_x, int start_y, int length_x, int length_y);

    /**
     * @brief Copies data from a source array to a block of xy-planes.
     *
     * @param src       Source array containing the data to copy.
     * @param start_x   Starting x-index of the block.
     * @param start_y   Starting y-index of the block.
     * @param length_x  Number of x-planes to copy.
     * @param length_y  Number of y-planes to copy.
     */
    void copy_xy_block_from(double* src, int start_x, int start_y, int length_x, int length_y);

    /**
     * @brief Sets all elements of the field to a specified value.
     *
     * @param clear_value The value to set for all elements (default is 0.0).
     */
    void clear(double clear_value = 0.0);

    /**
     * @brief Swaps the contents of two field3 objects.
     *
     * @param lhs       The first field3 object.
     * @param rhs       The second field3 object.
     */
    friend void swap(field3& lhs, field3& rhs);

protected:
    int         nx = 0, ny = 0, nz = 0;
    int         size_n = 0;
    std::string name   = "Default";
};
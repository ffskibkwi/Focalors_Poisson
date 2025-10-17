#include "field3.h"

#include "base/pch.h"

/**
 * @brief Constructor with initialization parameters.
 *
 * Creates a field3 object and initializes it with the specified dimensions and name.
 *
 * @param in_nx     Number of grid points in the x-direction.
 * @param in_ny     Number of grid points in the y-direction.
 * @param in_nz     Number of grid points in the z-direction.
 * @param in_name   Name identifier for the field.
 */
field3::field3(int in_nx, int in_ny, int in_nz, const std::string& in_name) { init(in_nx, in_ny, in_nz, in_name); }

/**
 * @brief Destructor.
 *
 * Releases all allocated memory and resources.
 */
field3::~field3() { deinit(); }

/**
 * @brief Assignment operator.
 *
 * Copies the values from another field3 object to this one.
 *
 * @param rhs       The field3 object to copy from.
 * @return          A reference to this object.
 */
field3& field3::operator=(const field3& rhs)
{
    if (this != &rhs)
    {
        OPENMP_PARALLEL_FOR()
        for (int i = 0; i < get_nx(); i++)
        {
            for (int j = 0; j < get_ny(); j++)
            {
                for (int k = 0; k < get_nz(); k++)
                {
                    this->operator()(i, j, k) = rhs.operator()(i, j, k);
                }
            }
        }
    }
    return *this;
}

/**
 * @brief Initializes the field3 object.
 *
 * Allocates memory for the field data and sets the dimensions and name.
 *
 * @param in_nx     Number of grid points in the x-direction.
 * @param in_ny     Number of grid points in the y-direction.
 * @param in_nz     Number of grid points in the z-direction.
 * @param in_name   Name identifier for the field (optional).
 */
void field3::init(int in_nx, int in_ny, int in_nz, const std::string& in_name)
{
    nx   = in_nx;
    ny   = in_ny;
    nz   = in_nz;
    name = in_name;

    size_n = nx * ny * nz;
    if (value != nullptr)
    {
        delete[] value;
    }
    value = new double[size_n];

    clear(0.0);
}

/**
 * @brief Deinitializes the field3 object.
 *
 * Releases allocated memory and resets internal variables.
 */
void field3::deinit()
{
    if (value != nullptr)
    {
        delete[] value;
        value = nullptr;
    }
}

/**
 * @brief Adds a constant value to all elements in the field.
 *
 * @param t         The constant value to add.
 */
void field3::add(const double t)
{
    OPENMP_PARALLEL_FOR()
    for (int i = 0; i < nx; i++)
    {
        for (int j = 0; j < ny; j++)
        {
            for (int k = 0; k < nz; k++)
            {
                this->operator()(i, j, k) = this->operator()(i, j, k) + t;
            }
        }
    }
}

/**
 * @brief Transposes the field according to specified dimension permutation.
 *
 * Rearranges the data according to the permutation of dimensions.
 *
 * @param dst           The destination field for the transposed data.
 * @param permutation   Array specifying the new order of dimensions.
 */
void field3::transpose(field3& dst, const std::array<int, 3>& permutation)
{
    int nx = this->get_nx();
    int ny = this->get_ny();
    int nz = this->get_nz();

    std::array<int, 3> dims = {nx, ny, nz};
    dst.set_size(dims[permutation[0]], dims[permutation[1]], dims[permutation[2]]);

    OPENMP_PARALLEL_FOR()
    for (int i = 0; i < nx; i++)
    {
        for (int j = 0; j < ny; j++)
        {
            for (int k = 0; k < nz; k++)
            {
                std::array<int, 3> src_idx = {i, j, k};
                std::array<int, 3> dst_idx;
                for (int p = 0; p < 3; p++)
                {
                    dst_idx[p] = src_idx[permutation[p]];
                }
                dst(dst_idx[0], dst_idx[1], dst_idx[2]) = this->operator()(i, j, k);
            }
        }
    }
}

/**
 * @brief Applies an affine transformation to the field.
 *
 * Performs the operation: this = this + a*x + b for each element.
 *
 * @param a         The scaling factor for field x.
 * @param x         The field to scale and add.
 * @param b         The constant value to add.
 */
void field3::add_affine_transform(const double a, const field3& x, const double b)
{
    OPENMP_PARALLEL_FOR()
    for (int i = 0; i < nx; i++)
    {
        for (int j = 0; j < ny; j++)
        {
            for (int k = 0; k < nz; k++)
            {
                this->operator()(i, j, k) += a * x(i, j, k) + b;
            }
        }
    }
}

/**
 * @brief Computes the sum of all elements in the field.
 *
 * @return          The sum of all elements.
 */
double field3::sum()
{
    double sum = 0.;
    for (size_t i = 0; i < size_n; ++i)
        sum += value[i];
    return sum;
}

/**
 * @brief Computes the sum of elements in an xy-plane at a given z-index.
 *
 * @param k         The z-index of the plane.
 * @return          The sum of elements in the xy-plane.
 */
double field3::sum_at_xy_plane(int k)
{
    double sum = 0.;

    for (int i = 0; i < nx; i++)
    {
        for (int j = 0; j < ny; j++)
        {
            sum += value[i * ny * nz + j * nz + k];
        }
    }

    return sum;
}

/**
 * @brief Computes the sum of elements in an xz-plane at a given y-index.
 *
 * @param j         The y-index of the plane.
 * @return          The sum of elements in the xz-plane.
 */
double field3::sum_at_xz_plane(int j)
{
    double sum = 0.;

    for (int i = 0; i < nx; i++)
    {
        for (int k = 0; k < nz; k++)
        {
            sum += value[i * ny * nz + j * nz + k];
        }
    }

    return sum;
}

/**
 * @brief Computes the sum of elements in a yz-plane at a given x-index.
 *
 * @param i         The x-index of the plane.
 * @return          The sum of elements in the yz-plane.
 */
double field3::sum_at_yz_plane(int i)
{
    double sum = 0.;

    for (int j = 0; j < ny; j++)
    {
        for (int k = 0; k < nz; k++)
        {
            sum += value[i * ny * nz + j * nz + k];
        }
    }

    return sum;
}

/**
 * @brief Computes the mean value of all elements in the field.
 *
 * @return          The mean value of all elements.
 */
double field3::mean() { return sum() / (nx * ny * nz); }

/**
 * @brief Computes the mean value of elements in an xy-plane at a given z-index.
 *
 * @param k         The z-index of the plane.
 * @return          The mean value of elements in the xy-plane.
 */
double field3::mean_at_xy_plane(int k) { return sum_at_xy_plane(k) / (nx * ny); }

/**
 * @brief Computes the mean value of elements in an xz-plane at a given y-index.
 *
 * @param j         The y-index of the plane.
 * @return          The mean value of elements in the xz-plane.
 */
double field3::mean_at_xz_plane(int j) { return sum_at_xz_plane(j) / (nx * nz); }

/**
 * @brief Computes the mean value of elements in a yz-plane at a given x-index.
 *
 * @param i         The x-index of the plane.
 * @return          The mean value of elements in the yz-plane.
 */
double field3::mean_at_yz_plane(int i) { return sum_at_yz_plane(i) / (ny * nz); }

/**
 * @brief Accessor operator for modifying elements.
 *
 * Provides read/write access to a specific element in the field.
 *
 * @param i         The x-index.
 * @param j         The y-index.
 * @param k         The z-index.
 * @return          Reference to the element at the specified indices.
 */
double& field3::operator()(int i, int j, int k) { return value[i * ny * nz + j * nz + k]; }

/**
 * @brief Accessor operator for reading elements.
 *
 * Provides read-only access to a specific element in the field.
 *
 * @param i         The x-index.
 * @param j         The y-index.
 * @param k         The z-index.
 * @return          Value of the element at the specified indices.
 */
double field3::operator()(int i, int j, int k) const { return value[i * ny * nz + j * nz + k]; }

/**
 * @brief Sets the size of the field without reallocating memory.
 *
 * Changes the dimensions of the field if the total number of elements
 * does not exceed the current allocation.
 *
 * @param in_nx     New number of grid points in the x-direction.
 * @param in_ny     New number of grid points in the y-direction.
 * @param in_nz     New number of grid points in the z-direction.
 * @return          True if successful, false if reallocation would be needed.
 */
bool field3::set_size(int in_nx, int in_ny, int in_nz)
{
    if (in_nx * in_ny * in_nz <= size_n)
    {
        nx = in_nx;
        ny = in_ny;
        nz = in_nz;
        return true;
    }
    else
    {
        return false;
    }
}

/**
 * @brief Gets a pointer to a specific element in the field.
 *
 * @param i         The x-index.
 * @param j         The y-index.
 * @param k         The z-index.
 * @return          Pointer to the element at the specified indices.
 */
double* field3::get_ptr(int i, int j, int k) { return value + ny * nz * i + nz * j + k; }

/**
 * @brief Copies a slice of xz-plane data to a destination array.
 *
 * @param dest      Destination array for the copied data.
 * @param start     Starting y-index of the slice.
 * @param length    Number of y-planes to copy.
 */
void field3::copy_xz_slice_to(double* dest, int start, int length)
{
    for (int j = start; j < start + length; j++)
    {
        OPENMP_PARALLEL_FOR()
        for (int i = 0; i < get_nx(); i++)
        {
            for (int k = 0; k < get_nz(); k++)
            {
                dest[(j - start) * get_nx() * get_nz() + i * get_nz() + k] = this->operator()(i, j, k);
            }
        }
    }
}

/**
 * @brief Copies data from a source array to a slice of xz-planes.
 *
 * @param src       Source array containing the data to copy.
 * @param start     Starting y-index of the slice.
 * @param length    Number of y-planes to copy.
 */
void field3::copy_xz_slice_from(double* src, int start, int length)
{
    for (int j = start; j < start + length; j++)
    {
        OPENMP_PARALLEL_FOR()
        for (int i = 0; i < get_nx(); i++)
        {
            for (int k = 0; k < get_nz(); k++)
            {
                this->operator()(i, j, k) = src[(j - start) * get_nx() * get_nz() + i * get_nz() + k];
            }
        }
    }
}

/**
 * @brief Copies a block of xy-plane data to a destination array.
 *
 * @param dest      Destination array for the copied data.
 * @param start_x   Starting x-index of the block.
 * @param start_y   Starting y-index of the block.
 * @param length_x  Number of x-planes to copy.
 * @param length_y  Number of y-planes to copy.
 */
void field3::copy_xy_block_to(double* dest, int start_x, int start_y, int length_x, int length_y)
{
    for (int i = start_x; i < start_x + length_x; i++)
    {
        for (int j = start_y; j < start_y + length_y; j++)
        {
            for (int k = 0; k < get_nz(); k++)
            {
                dest[(i - start_x) * length_y * get_nz() + (j - start_y) * get_nz() + k] = this->operator()(i, j, k);
            }
        }
    }
}

/**
 * @brief Copies data from a source array to a block of xy-planes.
 *
 * @param src       Source array containing the data to copy.
 * @param start_x   Starting x-index of the block.
 * @param start_y   Starting y-index of the block.
 * @param length_x  Number of x-planes to copy.
 * @param length_y  Number of y-planes to copy.
 */
void field3::copy_xy_block_from(double* src, int start_x, int start_y, int length_x, int length_y)
{
    for (int i = start_x; i < start_x + length_x; i++)
    {
        for (int j = start_y; j < start_y + length_y; j++)
        {
            for (int k = 0; k < get_nz(); k++)
            {
                this->operator()(i, j, k) = src[(i - start_x) * length_y * get_nz() + (j - start_y) * get_nz() + k];
            }
        }
    }
}

/**
 * @brief Sets all elements of the field to a specified value.
 *
 * @param clear_value The value to set for all elements (default is 0.0).
 */
void field3::clear(double clear_value)
{
    OPENMP_PARALLEL_FOR()
    for (int i = 0; i < nx; i++)
    {
        for (int j = 0; j < ny; j++)
        {
            for (int k = 0; k < nz; k++)
            {
                value[i * ny * nz + j * nz + k] = clear_value;
            }
        }
    }
}

/**
 * @brief Swaps the contents of two field3 objects.
 *
 * @param lhs       The first field3 object.
 * @param rhs       The second field3 object.
 */
void swap(field3& lhs, field3& rhs)
{
    using std::swap;

    swap(lhs.value, rhs.value);
    swap(lhs.nx, rhs.nx);
    swap(lhs.ny, rhs.ny);
    swap(lhs.nz, rhs.nz);
    swap(lhs.size_n, rhs.size_n);
    swap(lhs.name, rhs.name);
}

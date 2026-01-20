#pragma once

#include <array>
#include <string>

class field3
{
public:
    double* value = nullptr;

    field3() {}
    /**
     * @brief Constructor with dimensions only.
     *
     * Creates a field3 object and initializes it with the specified dimensions.
     *
     * @param in_nx     Number of grid points in the x-direction.
     * @param in_ny     Number of grid points in the y-direction.
     * @param in_nz     Number of grid points in the z-direction.
     */
    field3(int in_nx, int in_ny, int in_nz);
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
    field3(int in_nx, int in_ny, int in_nz, const std::string& in_name);
    /**
     * @brief Constructor with name only.
     *
     * Creates a field3 object with only a name, without initializing memory.
     * Memory will be allocated when init() is called.
     *
     * @param in_name   Name identifier for the field.
     */
    field3(const std::string& in_name);
    /**
     * @brief Destructor.
     *
     * Releases all allocated memory and resources.
     */
    ~field3();
    /**
     * @brief Assignment operator.
     *
     * Copies the values from another field3 object to this one.
     *
     * @param rhs       The field3 object to copy from.
     * @return          A reference to this object.
     */
    field3& operator=(const field3& rhs);
    /**
     * @brief Move constructor.
     *
     * Efficiently transfers resources from another field3 object.
     *
     * @param rhs       The field3 object to move from.
     */
    field3(field3&& rhs) noexcept { swap(*this, rhs); }
    /**
     * @brief Move assignment operator.
     *
     * Efficiently transfers resources from another field3 object to this one.
     *
     * @param rhs       The field3 object to move from.
     * @return          A reference to this object.
     */
    field3& operator=(field3&& rhs) noexcept
    {
        if (this != &rhs)
        {
            swap(*this, rhs);
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
    void init(int in_nx, int in_ny, int in_nz, const std::string& in_name = "");
    /**
     * @brief Deinitializes the field3 object.
     *
     * Releases allocated memory and resets internal variables.
     */
    void deinit();

    /**
     * @brief Adds a constant value to all elements in the field.
     *
     * @param t         The constant value to add.
     */
    void add(const double t);
    /**
     * @brief Applies an affine transformation to the field.
     *
     * Performs the operation: this = this + a*x + b for each element.
     *
     * @param a         The scaling factor for field x.
     * @param x         The field to scale and add.
     * @param b         The constant value to add.
     */
    void add_affine_transform(const double a, const field3& x, const double b);
    /**
     * @brief Transposes the field according to specified dimension permutation.
     *
     * Rearranges the data according to the permutation of dimensions.
     *
     * @param dst           The destination field for the transposed data.
     * @param permutation   Array specifying the new order of dimensions.
     */
    void transpose(field3& dst, const std::array<int, 3>& permutation);

    double norm();

    /**
     * @brief Computes the sum of all elements in the field.
     *
     * @return          The sum of all elements.
     */
    double sum();
    /**
     * @brief Computes the sum of elements in an xy-plane at a given z-index.
     *
     * @param k         The z-index of the plane.
     * @return          The sum of elements in the xy-plane.
     */
    double sum_at_xy_plane(int k);
    /**
     * @brief Computes the sum of elements in an xz-plane at a given y-index.
     *
     * @param j         The y-index of the plane.
     * @return          The sum of elements in the xz-plane.
     */
    double sum_at_xz_plane(int j);
    /**
     * @brief Computes the sum of elements in a yz-plane at a given x-index.
     *
     * @param i         The x-index of the plane.
     * @return          The sum of elements in the yz-plane.
     */
    double sum_at_yz_plane(int i);
    /**
     * @brief Computes the mean value of all elements in the field.
     *
     * @return          The mean value of all elements.
     */
    double mean();
    /**
     * @brief Computes the mean value of elements in an xy-plane at a given z-index.
     *
     * @param k         The z-index of the plane.
     * @return          The mean value of elements in the xy-plane.
     */
    double mean_at_xy_plane(int k);
    /**
     * @brief Computes the mean value of elements in an xz-plane at a given y-index.
     *
     * @param j         The y-index of the plane.
     * @return          The mean value of elements in the xz-plane.
     */
    double mean_at_xz_plane(int j);
    /**
     * @brief Computes the mean value of elements in a yz-plane at a given x-index.
     *
     * @param i         The x-index of the plane.
     * @return          The mean value of elements in the yz-plane.
     */
    double mean_at_yz_plane(int i);

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
    int         nx, ny, nz;
    int         size_n;
    std::string name = "Default";
};
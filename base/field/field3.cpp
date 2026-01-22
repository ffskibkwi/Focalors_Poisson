#include "field3.h"

#include "base/pch.h"
#include "field_macro.h"

field3::field3(int in_nx, int in_ny, int in_nz) { init(in_nx, in_ny, in_nz); }

field3::field3(int in_nx, int in_ny, int in_nz, const std::string& in_name) { init(in_nx, in_ny, in_nz, in_name); }

field3::field3(const std::string& in_name) { name = in_name; }

field3::~field3() { deinit(); }

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

void field3::init(int in_nx, int in_ny, int in_nz, const std::string& in_name)
{
    ASSERT_FIELD3_POSITIVE(in_nx, in_ny, in_nz, name);

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

void field3::deinit()
{
    if (value != nullptr)
    {
        delete[] value;
        value = nullptr;
    }
}

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

field3 field3::operator+(const field3& rhs)
{
    field3 R(nx, ny, nz);

    OPENMP_PARALLEL_FOR()
    for (int i = 0; i < nx; i++)
        for (int j = 0; j < ny; j++)
            for (int k = 0; k < nz; k++)
                R(i, j, k) = this->operator()(i, j, k) + rhs.operator()(i, j, k);
    return R;
}

field3 field3::operator-(const field3& rhs)
{
    field3 R(nx, ny, nz);

    OPENMP_PARALLEL_FOR()
    for (int i = 0; i < nx; i++)
        for (int j = 0; j < ny; j++)
            for (int k = 0; k < nz; k++)
                R(i, j, k) = this->operator()(i, j, k) - rhs.operator()(i, j, k);
    return R;
}

field3 field3::operator*(const double a)
{
    field3 R(nx, ny, nz);

    OPENMP_PARALLEL_FOR()
    for (int i = 0; i < nx; i++)
        for (int j = 0; j < ny; j++)
            for (int k = 0; k < nz; k++)
                R(i, j, k) = this->operator()(i, j, k) * a;
    return R;
}

field3& field3::operator*=(const double a)
{
    OPENMP_PARALLEL_FOR()
    for (int i = 0; i < nx; i++)
        for (int j = 0; j < ny; j++)
            for (int k = 0; k < nz; k++)
                this->operator()(i, j, k) *= a;
    return *this;
}

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

double field3::norm()
{
    double sum = 0.;
    for (size_t i = 0; i < size_n; ++i)
        sum += value[i] * value[i];
    return std::sqrt(sum);
}

double field3::sum()
{
    double sum = 0.;
    for (size_t i = 0; i < size_n; ++i)
        sum += value[i];
    return sum;
}

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

double field3::mean() { return sum() / (nx * ny * nz); }

double field3::mean_at_xy_plane(int k) { return sum_at_xy_plane(k) / (nx * ny); }

double field3::mean_at_xz_plane(int j) { return sum_at_xz_plane(j) / (nx * nz); }

double field3::mean_at_yz_plane(int i) { return sum_at_yz_plane(i) / (ny * nz); }

void field3::left_bond_add(const double a, const field2& bound)
{
    for (int j = 0; j < ny; j++)
        for (int k = 0; k < nz; k++)
            this->operator()(0, j, k) += a * bound(j, k);
}

void field3::right_bond_add(const double a, const field2& bound)
{
    for (int j = 0; j < ny; j++)
        for (int k = 0; k < nz; k++)
            this->operator()(nx - 1, j, k) += a * bound(j, k);
}

void field3::front_bond_add(const double a, const field2& bound)
{
    for (int i = 0; i < nx; i++)
        for (int k = 0; k < nz; k++)
            this->operator()(i, 0, k) += a * bound(i, k);
}

void field3::back_bond_add(const double a, const field2& bound)
{
    for (int i = 0; i < nx; i++)
        for (int k = 0; k < nz; k++)
            this->operator()(i, ny - 1, k) += a * bound(i, k);
}

void field3::down_bond_add(const double a, const field2& bound)
{
    for (int i = 0; i < nx; i++)
        for (int j = 0; j < ny; j++)
            this->operator()(i, j, 0) += a * bound(i, j);
}

void field3::up_bond_add(const double a, const field2& bound)
{
    for (int i = 0; i < nx; i++)
        for (int j = 0; j < ny; j++)
            this->operator()(i, j, nz - 1) += a * bound(i, j);
}

void field3::bond_add(LocationType location, const double a, const field2& bound)
{
    switch (location)
    {
        case LocationType::Left:
            left_bond_add(a, bound);
            break;
        case LocationType::Right:
            right_bond_add(a, bound);
            break;
        case LocationType::Front:
            front_bond_add(a, bound);
            break;
        case LocationType::Back:
            back_bond_add(a, bound);
            break;
        case LocationType::Down:
            down_bond_add(a, bound);
            break;
        case LocationType::Up:
            up_bond_add(a, bound);
            break;
        default:
            throw std::invalid_argument("Invalid location type");
    }
}

void field3::left_bond_add(const double a, const field3& neighbour)
{
    int neighbour_nx = neighbour.get_nx();
    for (int j = 0; j < ny; j++)
        for (int k = 0; k < nz; k++)
            this->operator()(0, j, k) += a * neighbour(neighbour_nx - 1, j, k);
}

void field3::right_bond_add(const double a, const field3& neighbour)
{
    for (int j = 0; j < ny; j++)
        for (int k = 0; k < nz; k++)
            this->operator()(nx - 1, j, k) += a * neighbour(0, j, k);
}

void field3::front_bond_add(const double a, const field3& neighbour)
{
    int neighbour_ny = neighbour.get_ny();
    for (int i = 0; i < nx; i++)
        for (int k = 0; k < nz; k++)
            this->operator()(i, 0, k) += a * neighbour(i, neighbour_ny - 1, k);
}

void field3::back_bond_add(const double a, const field3& neighbour)
{
    for (int i = 0; i < nx; i++)
        for (int k = 0; k < nz; k++)
            this->operator()(i, ny - 1, k) += a * neighbour(i, 0, k);
}

void field3::down_bond_add(const double a, const field3& neighbour)
{
    int neighbour_nz = neighbour.get_nz();
    for (int i = 0; i < nx; i++)
        for (int j = 0; j < ny; j++)
            this->operator()(i, j, 0) += a * neighbour(i, j, neighbour_nz - 1);
}

void field3::up_bond_add(const double a, const field3& neighbour)
{
    for (int i = 0; i < nx; i++)
        for (int j = 0; j < ny; j++)
            this->operator()(i, j, nz - 1) += a * neighbour(i, j, 0);
}

void field3::bond_add(LocationType location, const double a, const field3& neighbour)
{
    switch (location)
    {
        case LocationType::Left:
            left_bond_add(a, neighbour);
            break;
        case LocationType::Right:
            right_bond_add(a, neighbour);
            break;
        case LocationType::Front:
            front_bond_add(a, neighbour);
            break;
        case LocationType::Back:
            back_bond_add(a, neighbour);
            break;
        case LocationType::Down:
            down_bond_add(a, neighbour);
            break;
        case LocationType::Up:
            up_bond_add(a, neighbour);
            break;
        default:
            throw std::invalid_argument("Invalid location type");
    }
}

double& field3::operator()(int i, int j, int k)
{
    ASSERT_FIELD3_BOUNDS(i, j, k, nx, ny, nz, name);
    return value[i * ny * nz + j * nz + k];
}

double field3::operator()(int i, int j, int k) const
{
    ASSERT_FIELD3_BOUNDS(i, j, k, nx, ny, nz, name);
    return value[i * ny * nz + j * nz + k];
}

double* field3::get_ptr(int i, int j, int k) const
{
    ASSERT_FIELD3_BOUNDS(i, j, k, nx, ny, nz, name);
    return value + ny * nz * i + nz * j + k;
}

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

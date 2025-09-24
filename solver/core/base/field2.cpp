#include "field2.h"

#include "pch.h"

field2::field2(int in_nx, int in_ny, const std::string& in_name) { init(in_nx, in_ny, in_name); }

field2::~field2() { deinit(); }

field2::field2(const field2& rhs) noexcept
{
    if (value != nullptr)
    {
        delete[] value;
        value = nullptr;
    }

    nx     = rhs.nx;
    ny     = rhs.ny;
    size_n = rhs.size_n;
    value  = new double[size_n];
    for (int i = 0; i < size_n; i++)
        value[i] = rhs.value[i];
}

field2& field2::operator=(const field2& rhs) noexcept
{
    if (this == &rhs)
        return *this;

    if (value != nullptr)
    {
        delete[] value;
        value = nullptr;
    }

    nx     = rhs.nx;
    ny     = rhs.ny;
    size_n = rhs.size_n;
    value  = new double[size_n];
    for (int i = 0; i < size_n; i++)
        value[i] = rhs.value[i];

    return *this;
}

void field2::init(int in_nx, int in_ny, const std::string& in_name)
{
    nx   = in_nx;
    ny   = in_ny;
    name = in_name;

    size_n = nx * ny;
    if (value != nullptr)
    {
        delete[] value;
    }
    value = new double[size_n];

    clear(0.0);
}

void field2::deinit()
{
    if (value != nullptr)
    {
        delete[] value;
        value = nullptr;
    }
}

field2 field2::operator+(const field2& rhs)
{
    field2 R(nx, ny);

    OPENMP_PARALLEL_FOR()
    for (int i = 0; i < nx; i++)
    {
        for (int j = 0; j < ny; j++)
        {
            R(i, j) = this->operator()(i, j) + rhs.operator()(i, j);
        }
    }
    return R;
}

field2 field2::operator-(const field2& rhs)
{
    field2 R(nx, ny);

    OPENMP_PARALLEL_FOR()
    for (int i = 0; i < nx; i++)
    {
        for (int j = 0; j < ny; j++)
        {
            R(i, j) = this->operator()(i, j) - rhs.operator()(i, j);
        }
    }
    return R;
}

field2 field2::operator*(const double k)
{
    field2 R(nx, ny);

    OPENMP_PARALLEL_FOR()
    for (int i = 0; i < nx; i++)
    {
        for (int j = 0; j < ny; j++)
        {
            R(i, j) = this->operator()(i, j) * k;
        }
    }
    return R;
}

field2 field2::operator*(const field2& rhs)
{
    int m = this->nx;
    int n = this->ny;
    int p = rhs.ny;

    field2 R(m, p);

    OPENMP_PARALLEL_FOR()
    for (int i = 0; i < m; i++)
    {
        for (int j = 0; j < p; j++)
        {
            R(i, j) = 0.0;
            for (int k = 0; k < n; k++)
            {
                R(i, j) += this->operator()(i, k) * rhs(k, j);
            }
        }
    }
    return R;
}

void field2::add(const double t)
{
    OPENMP_PARALLEL_FOR()
    for (int i = 0; i < nx; i++)
    {
        for (int j = 0; j < ny; j++)
        {
            this->operator()(i, j) = this->operator()(i, j) + t;
        }
    }
}

void field2::add_affine_transform(const double a, const field2& x, const double b)
{
    OPENMP_PARALLEL_FOR()
    for (int i = 0; i < nx; i++)
    {
        for (int j = 0; j < ny; j++)
        {
            this->operator()(i, j) += a * x(i, j) + b;
        }
    }
}

double field2::dot(const field2& A)
{
    double sum = 0.;
    for (int i = 0; i < size_n; i++)
        sum += value[i] * A.value[i];
    return sum;
}

double field2::norm()
{
    double sum = 0.;
    for (size_t i = 0; i < size_n; ++i)
        sum += value[i] * value[i];
    return std::sqrt(sum);
}

double field2::sum()
{
    double sum = 0.;
    for (size_t i = 0; i < size_n; ++i)
        sum += value[i];
    return sum;
}

double field2::squared_sum()
{
    double sum = 0.;
    for (size_t i = 0; i < size_n; ++i)
        sum += value[i] * value[i];
    return sum;
}

double field2::mean_at_x_axis(int i)
{
    double sum = 0;
    for (int j = 0; j < ny; j++)
    {
        sum += value[i * ny + j];
    }
    return sum / ny;
}

double field2::mean_at_y_axis(int j)
{
    double sum = 0;
    for (int i = 0; i < nx; i++)
    {
        sum += value[i * ny + j];
    }
    return sum / nx;
}

void field2::left_bond_add(const double k, const field2& left_field2)
{
    int left_field2_nx = left_field2.get_nx();
    for (int j = 0; j < ny; j++)
        value[j] += k * left_field2(left_field2_nx - 1, j);
}

void field2::right_bond_add(const double k, const field2& right_field2)
{
    // int right_field2_nx = right_field2.get_nx();
    for (int j = 0; j < ny; j++)
        value[(nx - 1) * ny + j] += k * right_field2(0, j);
}

void field2::up_bond_add(const double k, const field2& up_field2)
{
    // int up_field2_ny = up_field2.get_ny();
    for (int i = 0; i < nx; i++)
        value[i * ny + ny - 1] += k * up_field2(i, 0);
}

void field2::down_bond_add(const double k, const field2& down_field2)
{
    int down_field2_ny = down_field2.get_ny();
    for (int i = 0; i < nx; i++)
        value[i * ny] += k * down_field2(i, down_field2_ny - 1);
}

double& field2::operator()(int i, int j) { return value[i * ny + j]; }

double field2::operator()(int i, int j) const { return value[i * ny + j]; }

bool field2::set_size(int in_nx, int in_ny)
{
    if (in_nx * in_ny <= size_n)
    {
        nx = in_nx;
        ny = in_ny;
        return true;
    }
    else
    {
        return false;
    }
}

double* field2::get_ptr(int i, unsigned int j) { return value + ny * i + j; }

double* field2::get_ptr(int i) { return value + ny * i; }

void field2::clear(double clear_value)
{
    OPENMP_PARALLEL_FOR()
    for (int i = 0; i < nx; i++)
    {
        for (int j = 0; j < ny; j++)
        {
            value[i * ny + j] = 0.0;
        }
    }
}

void swap(field2& lhs, field2& rhs)
{
    using std::swap;

    swap(lhs.value, rhs.value);
    swap(lhs.nx, rhs.nx);
    swap(lhs.ny, rhs.ny);
    swap(lhs.size_n, rhs.size_n);
    swap(lhs.name, rhs.name);
}

/**
 * @brief Transposes the field.
 *
 * Rearranges the data.
 *
 * @param dst           The destination field for the transposed data.
 */
void field2::transpose(field2& dst)
{
    int nx = this->get_nx();
    int ny = this->get_ny();

    dst.set_size(ny, nx);

    OPENMP_PARALLEL_FOR()
    for (int i = 0; i < nx; i++)
    {
        for (int j = 0; j < ny; j++)
        {
            dst(j, i) = this->operator()(i, j);
        }
    }
}
#include "field2.h"

#include "base/pch.h"
#include "field_macro.h"

field2::field2(int in_nx, int in_ny, const std::string& in_name) { init(in_nx, in_ny, in_name); }

field2::field2(const std::string& in_name) { name = in_name; }

field2::field2(int in_nx, int in_ny) { init(in_nx, in_ny); }

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
    init(in_nx, in_ny);
    name = in_name;
}

void field2::init(int in_nx, int in_ny)
{
    ASSERT_FIELD2_POSITIVE(in_nx, in_ny, name);

    nx = in_nx;
    ny = in_ny;

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
        for (int j = 0; j < ny; j++)
            R(i, j) = this->operator()(i, j) + rhs.operator()(i, j);
    return R;
}

field2 field2::operator-(const field2& rhs)
{
    field2 R(nx, ny);

    OPENMP_PARALLEL_FOR()
    for (int i = 0; i < nx; i++)
        for (int j = 0; j < ny; j++)
            R(i, j) = this->operator()(i, j) - rhs.operator()(i, j);
    return R;
}

field2 field2::operator*(const double a)
{
    field2 R(nx, ny, this->get_name());

    OPENMP_PARALLEL_FOR()
    for (int i = 0; i < nx; i++)
        for (int j = 0; j < ny; j++)
            R(i, j) = this->operator()(i, j) * a;
    return R;
}

field2 field2::operator*(const field2& rhs)
{
    int m = this->nx;
    int n = this->ny;
    int p = rhs.ny;

    field2 R(m, p, this->get_name() + '*' + rhs.get_name());

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

field2& field2::operator*=(const double a)
{
    OPENMP_PARALLEL_FOR()
    for (int i = 0; i < nx; i++)
        for (int j = 0; j < ny; j++)
            this->operator()(i, j) *= a;
    return *this;
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

void field2::left_bond_add(const double a, double* bound)
{
    for (int j = 0; j < ny; j++)
        this->operator()(0, j) += a * bound[j];
}

void field2::right_bond_add(const double a, double* bound)
{
    for (int j = 0; j < ny; j++)
        this->operator()(nx - 1, j) += a * bound[j];
}

void field2::down_bond_add(const double a, double* bound)
{
    for (int i = 0; i < nx; i++)
        this->operator()(i, 0) += a * bound[i];
}

void field2::up_bond_add(const double a, double* bound)
{
    for (int i = 0; i < nx; i++)
        this->operator()(i, ny - 1) += a * bound[i];
}

void field2::bond_add(LocationType location, const double a, double* bound)
{
    switch (location)
    {
        case LocationType::Left:
            left_bond_add(a, bound);
            break;
        case LocationType::Right:
            right_bond_add(a, bound);
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

void field2::left_bond_add(const double a, const field2& neighbour)
{
    int neighbour_nx = neighbour.get_nx();
    for (int j = 0; j < ny; j++)
        this->operator()(0, j) += a * neighbour(neighbour_nx - 1, j);
}

void field2::right_bond_add(const double a, const field2& neighbour)
{
    for (int j = 0; j < ny; j++)
        this->operator()(nx - 1, j) += a * neighbour(0, j);
}

void field2::down_bond_add(const double a, const field2& neighbour)
{
    int neighbour_ny = neighbour.get_ny();
    for (int i = 0; i < nx; i++)
        this->operator()(i, 0) += a * neighbour(i, neighbour_ny - 1);
}

void field2::up_bond_add(const double a, const field2& neighbour)
{
    for (int i = 0; i < nx; i++)
        this->operator()(i, ny - 1) += a * neighbour(i, 0);
}

void field2::bond_add(LocationType location, const double a, const field2& neighbour)
{
    std::cout << "before bond add " << location << std::endl;
    std::cout << "neighbour " << neighbour.get_name() << std::endl;
    neighbour.print();

    switch (location)
    {
        case LocationType::Left:
            left_bond_add(a, neighbour);
            break;
        case LocationType::Right:
            right_bond_add(a, neighbour);
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

    std::cout << "after bond add " << location << std::endl;
    std::cout << "neighbour " << neighbour.get_name() << std::endl;
    neighbour.print();
}

double& field2::operator()(int i, int j)
{
    ASSERT_FIELD2_BOUNDS(i, j, nx, ny, name);
    return value[i * ny + j];
}

double field2::operator()(int i, int j) const
{
    ASSERT_FIELD2_BOUNDS(i, j, nx, ny, name);
    return value[i * ny + j];
}

double* field2::get_ptr(int i, int j) const
{
    ASSERT_FIELD2_BOUNDS(i, j, nx, ny, name);
    return value + ny * i + j;
}

double* field2::get_ptr(int i) const
{
    ASSERT_BOUNDS(i, size_n, name);
    return value + ny * i;
}

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

void field2::clear(double clear_value)
{
    OPENMP_PARALLEL_FOR()
    for (int i = 0; i < nx; i++)
    {
        for (int j = 0; j < ny; j++)
        {
            value[i * ny + j] = clear_value;
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

void field2::print() const
{
    for (int i = 0; i < nx; i++)
    {
        for (int j = 0; j < ny; j++)
        {
            std::cout << this->operator()(i, j);
            if (j == ny - 1)
            {
                std::cout << '\n';
            }
            else
            {
                std::cout << ',';
            }
        }
    }
}
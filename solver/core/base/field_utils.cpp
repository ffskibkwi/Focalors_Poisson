#include "field_utils.hpp"

#include <array>

namespace FieldUtils
{
    void transpose(const field2& src, field2& dst)
    {
        int nx = src.get_nx();
        int ny = src.get_ny();

        dst.set_size(ny, nx);

        OPENMP_PARALLEL_FOR()
        for (int i = 0; i < nx; i++)
        {
            for (int j = 0; j < ny; j++)
            {
                dst(j, i) = src(i, j);
            }
        }
    }

    void transpose(const field3& src, field3& dst, const std::array<int, 3>& permutation)
    {
        int nx = src.get_nx();
        int ny = src.get_ny();
        int nz = src.get_nz();

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
                    dst(dst_idx[0], dst_idx[1], dst_idx[2]) = src(i, j, k);
                }
            }
        }
    }
} // namespace FieldUtils

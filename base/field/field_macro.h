#pragma once
#include <cassert>
#include <cstdio>
#include <string>

#ifndef NDEBUG

inline void assert_msg(bool cond, const std::string& msg)
{
    if (!cond)
    {
        std::fprintf(stderr, "Assertion failed: %s\n", msg.c_str());
        assert(cond);
    }
}

#    define ASSERT_BOUNDS(i, nx, name) assert_msg((i) >= 0 && (i) < (nx), name + ": i out of bounds");

// ---- 2D field checks ----
#    define ASSERT_FIELD2_POSITIVE(i, j, name)          \
        assert_msg((i) >= 0, name + ": i >= 0 failed"); \
        assert_msg((j) >= 0, name + ": j >= 0 failed")

#    define ASSERT_FIELD2_BOUNDS(i, j, nx, ny, name)                    \
        assert_msg((i) >= 0 && (i) < (nx), name + ": i out of bounds"); \
        assert_msg((j) >= 0 && (j) < (ny), name + ": j out of bounds")

// ---- 3D field checks ----
#    define ASSERT_FIELD3_POSITIVE(i, j, k, name)       \
        assert_msg((i) >= 0, name + ": i >= 0 failed"); \
        assert_msg((j) >= 0, name + ": j >= 0 failed"); \
        assert_msg((k) >= 0, name + ": k >= 0 failed")

#    define ASSERT_FIELD3_BOUNDS(i, j, k, nx, ny, nz, name)             \
        assert_msg((i) >= 0 && (i) < (nx), name + ": i out of bounds"); \
        assert_msg((j) >= 0 && (j) < (ny), name + ": j out of bounds"); \
        assert_msg((k) >= 0 && (k) < (nz), name + ": k out of bounds")

#else // NDEBUG defined → disable all asserts

#    define ASSERT_BOUNDS(i, nx, name)                      ((void)0)
#    define ASSERT_FIELD2_POSITIVE(i, j, name)              ((void)0)
#    define ASSERT_FIELD2_BOUNDS(i, j, nx, ny, name)        ((void)0)
#    define ASSERT_FIELD3_POSITIVE(i, j, k, name)           ((void)0)
#    define ASSERT_FIELD3_BOUNDS(i, j, k, nx, ny, nz, name) ((void)0)

#endif

#pragma once

#ifdef HYBRID
#    define STRINGIFY(x) #x
#    define STR(x)       STRINGIFY(x)
#    define OMP_WITH_CLAUSE(...) omp parallel for schedule(dynamic) __VA_ARGS__
#    define OPENMP_PARALLEL_FOR(...) _Pragma(STR(OMP_WITH_CLAUSE(__VA_ARGS__)))
#else
#    define OPENMP_PARALLEL_FOR(...) // No-op
#endif
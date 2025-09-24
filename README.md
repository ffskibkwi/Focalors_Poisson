## Env

Add fftw path in `~/.bashrc`

```shell
export FFTW_ROOT=/path/to/fftw3
```

## Build

```shell
cmake -S . -B build -DCMAKE_CXX_COMPILER=mpicxx -DCMAKE_BUILD_TYPE=Hybrid
cmake --build build --parallel $(nproc)
```
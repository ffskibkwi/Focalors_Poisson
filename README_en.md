<div align="center">
  <img src="Focalors_logo.png" alt="Focalors Logo" width="250"/>
</div>
<p align="center">
  # Focalors_Poisson
  <br/>
  High-Performance Poisson Solver Based on FFT-based Domain Decomposition Algorithm
</p>
<p align="center">
  [![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
  [![Language](https://img.img.shields.io/badge/Language-C%2B%2B-blue.svg)](https://isocpp.org/)
  [![CMake](https://img.img.shields.io/badge/Build-CMake-green.svg)](https://cmake.org/)
</p>
<p align="center">
  [中文](./README.md)&nbsp;&nbsp;|&nbsp;&nbsp;English&nbsp;&nbsp;|&nbsp;&nbsp;[日本語](README_ja.md)
</p>

**Focalors_Poisson** is a high-performance Poisson equation solver designed for composite geometries. It is based on the Fast Fourier Transform (FFT) and Domain Decomposition Method, aiming to deliver extreme computational speed.

This project serves as the Poisson module developed for the (work-in-progress) incompressible Navier-Stokes solver [Focalors](), but it can also operate as a fully functional, stand-alone Poisson equation solver.

---

## Core Philosophy

Our first principle of development is: **Speed** (potentially the world's fastest Poisson equation solver for complex composite geometries).

*   **⚡️ Low Algorithmic Complexity**: The complexity of the entire solution process is controlled to $O(N \log N)$, which is comparable to that of the FFT algorithm itself.
*   **🚀 Highly Parallel**: The entire computational workflow of the solver fully supports parallel computing. (The MPI version is currently undergoing testing and optimization).
*   **💻 Rapid Development**: Highly encapsulated data structures and semantically defined solution procedures allow users to avoid tedious, repetitive configurations and focus on the solution itself.
*   **🔧 Easy to Use**: Only a standard C++ compiler and the FFTW library are required to effortlessly compile and run this project, eliminating the need for complex dependency configurations.

## Algorithm and Numerical Validation

The core algorithms employed by this solver and their detailed numerical validations can be found in:

*   **arXiv**: [https://arxiv.org/abs/2509.23180](https://arxiv.org/abs/2509.23180)
*   **DOI**: [https://doi.org/10.48550/arXiv.2509.23180](https://doi.org/10.48550/arXiv.2509.23180)

## Getting Started

Please follow the steps below to compile and run Focalors in your local environment.

### 1. Prerequisites

All development and testing are currently performed on Linux, so **it is recommended to use this solver on a Linux system**, ensuring that the system has the following software/packages:

*   **C++ Compiler**: Such as GCC (g++) or Intel C++ Compiler.
*   **FFTW3**: Currently the most well-known high-performance Fast Fourier Transform library.
    *   Please download the source code and compile it from the [FFTW official website](http://www.fftw.org/download.html).
*   **CMake**

### 2. Configuration

First, clone this repository to your local machine:
```shell
git clone https://github.com/ffskibkwi/Focalors_Poisson.git
cd Focalors_Poisson
```

Next, configure the environment variable to specify the installation path of the FFTW library:
```shell
export FFTW_ROOT=/path/to/fftw3 # Replace /path/to/fftw3 with the actual installation path of FFTW.
```

### 3. Compilation

(You can complete your personal code writing before this step)

```shell
cmake -S . -B build -DCMAKE_CXX_COMPILER=mpicxx -DCMAKE_BUILD_TYPE=Release

# 2. Compile in parallel using multiple cores
cmake --build build --parallel $(nproc) # Parallel compilation is recommended; replace $(nproc) with the number of cores on your machine.
```

### 4. Running

The compiled executable files are located in the `build/bin/` directory. Simply run them to start solving.

## How to Use

To use Focalors_Poisson to solve a specific problem, you need to create a `case` file – a C++ source file used to define the geometry, boundary conditions, and solution parameters.

**It is highly recommended to refer to `case/validation_cases/cross_new.cpp` as an example to start writing.**

**Basic Steps**:
1.  Create a new `.cpp` file in the `case/` directory.
2.  Include Focalors' core header files.
3.  Define your computational domain and composite geometry.
4.  Set boundary conditions and the right-hand side (source term).
5.  Initialize and invoke the solver.
6.  Output the results to a file.
7.  (By default, the `CMakeLists.txt` in the root directory already includes all files under `case/`, so compilation can proceed directly.)

(Detailed development documentation is still actively being produced.)

## Contributing

We welcome contributions in any form! If you find a bug, have feature suggestions, or wish to improve the code, please feel free to:
1.  Submit an [Issue](https://github.com/your-username/Focalors/issues).
2.  Fork this repository and open a Pull Request.

## Developers

The core algorithm design and code implementation of this solver were exclusively completed by the following two individuals:

*   **Zichao Jiang** - [jzc_focal@hotmail.com](mailto:jzc_focal@hotmail.com)
*   **Jiacheng Lian** - [blizzzzzzzzzzzzard@gmail.com](mailto:blizzzzzzzzzzzzard@gmail.com)

## Why "Focalors"?

Yes, it refers to the Hydro Archon from "Genshin Impact".

## License

This project is licensed under the MIT License. Please see the [LICENSE](LICENSE) file for details.
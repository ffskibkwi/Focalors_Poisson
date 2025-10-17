<p align="center">
  <img src="Focalors_logo.png" alt="Focalors Logo" width="250"/>
</p>
<div align="center">

  # Focalors_Poisson

  基于FFT区域分解算法的高性能Poisson求解器
</div>

<div align="center">

  [![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
  [![Language](https://img.shields.io/badge/Language-C%2B%2B-blue.svg)](https://isocpp.org/)
  [![CMake](https://img.shields.io/badge/Build-CMake-green.svg)](https://cmake.org/)

  中文&nbsp;&nbsp;|&nbsp;&nbsp;[English](./README_en.md)&nbsp;&nbsp;|&nbsp;&nbsp;[日本語](./README_ja.md)

</div>

---

**Focalors_Poisson** 是为拼接几何问题设计的高性能Poisson方程求解器。它基于快速傅里叶变换（FFT）和区域分解算法（Domain Decomposition Method），旨在提供最快的计算速度。

本项目是为（开发中的）不可压缩 Navier-Stokes 求解器 [Focalors](https://github.com/ffskibkwi/Focalors) 开发的 Poisson 模块，但它也可以作为一个功能完整、完全独立的Poisson方程求解器运行。

## 核心理念

开发的第一原则：**快**（也许是世界上最快的、适用于复杂拼接几何的Poisson方程求解器）

*   **⚡️ 低算法复杂度**: 整个求解流程的复杂度控制在了与FFT算法相当的$O\left(N \log N\right)$复杂度
*   **🚀 高度并行**: 求解器的整个计算流程完全支持并行计算（MPI版本目前正在测试与优化中）
*   **💻 快速开发**: 高度封装的数据结构和语义化的求解流程定义，让使用者可以避免繁琐的重复配置，专注于求解本身
*   **🔧 便捷使用**: 仅需标准的C++编译器和FFTW库，即可轻松编译和运行本项目，无需复杂的依赖配置

## 算法与数值验证

本求解器所采用的核心算法及其详细的数值验证参见：

*   **arXiv**: [https://arxiv.org/abs/2509.23180](https://arxiv.org/abs/2509.23180)
*   **DOI**: [https://doi.org/10.48550/arXiv.2509.23180](https://doi.org/10.48550/arXiv.2509.23180)

## 快速入门

请按照以下步骤在本地环境中编译并运行 Focalors

### 1. 环境准备

目前所有的开发和测试均在Linux上完成，所以**推荐在Liunx系统上使用该求解器**，并确保系统已具备以下软件/包：

*   **C++ 编译器**: GCC (g++) 或 Intel C++ Complier 等
*   **FFTW3**: 目前最知名的高效快速傅里叶变换库
    *   请从 [FFTW 官网](http://www.fftw.org/download.html) 下载源码编译安装
*   **CMake**

### 2. 配置

首先，克隆本仓库到本地：
```shell
git clone https://github.com/ffskibkwi/Focalors_Poisson.git
cd Focalors_Poisson
```

接下来，配置环境变量，指定 FFTW 库的安装路径：
```shell
export FFTW_ROOT=/path/to/fftw3 # 将 /path/to/fftw3 替换为FFTW的实际安装路径
```

### 3. 编译
（在本步骤前可以在完成个人代码的编写）

```shell
cmake -S . -B build -DCMAKE_CXX_COMPILER=mpicxx -DCMAKE_BUILD_TYPE=Release
cmake --build build --parallel $(nproc) #推荐使用并行编译，将$(nproc)替换为机器的核心数
```

### 4. 运行

编译生成的可执行文件位于 `build/bin/` 目录下，直接运行即可开始求解

## 如何使用

要使用 Focalors_Poisson 求解特定的问题，您需要创建一个 `case` 文件 -- 一个用于定义几何、边界条件和求解参数的 C++ 源文件。

**强烈建议参考 `case/validation_cases/cross_new.cpp` 作为示例来开始编写**

**基本步骤**:
1.  在 `case/` 目录下创建一个新的 `.cpp` 文件
2.  包含 Focalors 的核心头文件
3.  定义拼接几何`Geometry`与边界条件
4.  在构建的`Variable`中构建右端项
5.  初始化并调用求解器
6.  将结果输出到文件
7.  （默认情况下，根目录的 `CMakeLists.txt` 中已经包含了`case/`目录下的所有文件，所以可以直接进行编译）

（详细的开发文档仍在锐意制作中）

## 贡献

我们欢迎任何形式的贡献！如果您发现 Bug、有功能建议或希望改进代码，请随时：
1.  提交一个 [Issue](https://github.com/your-username/Focalors/issues)。
2.  Fork 本仓库并发起一个 Pull Request。

## 开发者

本求解器的核心算法设计与代码实现由且仅有以下两人完成：

*   **Zichao Jiang** - [jzc_focal@hotmail.com](mailto:jzc_focal@hotmail.com)
*   **Jiacheng Lian** - [blizzzzzzzzzzzzard@gmail.com](mailto:blizzzzzzzzzzzzard@gmail.com)

## 为什么叫Focalors？
是的，就是《原神》中水之国度的神明~

## 许可证

本项目采用 MIT 许可证。详情请见 [LICENSE](LICENSE) 文件。
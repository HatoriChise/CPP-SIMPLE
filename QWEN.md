# CFD SIMPLE Solver - 项目上下文

## 项目概述

这是一个基于SIMPLE（Semi-Implicit Method for Pressure-Linked Equations）算法的二维CFD求解器，使用有限体积法离散化Navier-Stokes方程。项目采用现代C++17开发，使用结构化网格处理二维流体力学问题，如顶盖驱动流（lid-driven cavity flow）。

## 技术栈

- **语言**: C++17
- **构建系统**: CMake 3.14+
- **核心依赖**:
  - `Eigen3`: 线性代数运算（矩阵求解器）
  - `Boost`: `multi_array`（多维数组）、`math`（数学常量）
  - `fmt`: 格式化输出

## 项目结构

```
CPP-SIMPLE/
├── src/
│   ├── config/              # 配置参数
│   │   └── simulationConfig.hpp  # 全局仿真参数
│   ├── grid/                # 网格模块
│   │   ├── structuredMesh.hpp
│   │   └── structuredMesh.cpp
│   ├── field/               # 场数据模块
│   │   ├── scalarField.hpp/cpp    # 标量场
│   │   ├── vectorField.hpp/cpp    # 矢量场（速度）
│   │   ├── boundaryField.hpp/cpp  # 边界条件场
│   │   └── fluidProperties.hpp/cpp # 流体物性场
│   ├── equation/            # 方程离散模块
│   │   ├── scalarEquation.hpp/cpp # 标量输运方程
│   ├── main.cpp             # 程序入口
│   └── test.cpp/hpp         # 测试和示例
├── utils/                   # 工具函数
│   └── formatter4eigen.h    # Eigen矩阵的fmt格式化支持
├── data/                    # 输出数据目录
├── build/                   # 构建目录（已忽略）
└── CMakeLists.txt           # 构建配置
```

## 构建和运行

### 构建命令

```bash
cd build
cmake ..
make
```

### 运行程序

```bash
./test_deps
```

### 依赖安装（vcpkg）

项目使用vcpkg管理依赖，需要安装：

- `fmt`
- `eigen3`
- `boost`

## 编码规范

### 文件组织

- 头文件使用 `#pragma once` 防护
- 类声明在 `.hpp` 文件，实现在 `.cpp` 文件
- 包含路径相对于项目根目录：`#include "src/config/simulationConfig.hpp"`

### 命名约定

- **类名**: PascalCase（如 `ScalarField`, `StructuredMesh`）
- **函数/方法**: camelCase（如 `createVolumeMesh`, `addDiffusionTerm`）
- **成员变量**: 带下划线后缀（如 `ncx_`, `data_`）
- **常量**: camelCase 或 UPPER_CASE（如 `ncx`, `Lx`, `Re`）
- **枚举**: UPPER_CASE（如 `FACE_POSITION`, `BOUNDARY_TYPE_U`）

### 数据布局

- 使用 `boost::multi_array<float, 2>` 存储二维场数据
- 索引映射：`operator()(i, j)` → `data_[j][i]`（x方向内存连续）
- 全局坐标：x从左到右（0到ncx-1），y从下到上（0到ncy-1）

### 代码风格

- 使用 `fmt` 库进行格式化输出
- 使用 `constexpr` 定义编译时常量
- 配置参数集中在 `simulationConfig.hpp` 中

## 核心模块

### 配置模块 (`src/config/simulationConfig.hpp`)

包含所有全局配置参数：

- 几何参数：`Lx`, `Ly`, `ncx`, `ncy`
- 物理参数：`Re`, `density`, `mu_fluid`, `k_fluid`, `cp_fluid`
- 边界条件：`boundaryInfo[4]`
- 求解控制：`max_iterations`, `continuity_tolerance`

### 场数据模块 (`src/field/`)

- `ScalarField`: 标量场（温度、压力等）
- `VectorField`: 矢量场（速度），包含u和v分量
- `BoundaryField`: 边界条件存储和访问
- `FluidPropertyField`: 流体物性（密度、粘度、导热率、比热）

### 方程离散模块 (`src/equation/`)

`ScalarEquation`: 标量输运方程离散化，包含：

- 系数结构 `COEF`: `aE`, `aW`, `aN`, `aS`, `aP`, `bsrc`
- 方法：`resetCoefficients()`, `addDiffusionTerm()`, `addConvectionTerm()` 等

## 开发注意事项

1. **边界处理**: `ScalarField`和`VectorField`不处理边界，边界条件由`BoundaryField`管理
2. **物性更新**: 当前使用常物性，`updateFluidProperties()`待实现变物性逻辑
3. **方程离散**: `ScalarEquation`的部分方法待实现
4. **内存管理**: `boost::multi_array`自动管理内存

## SIMPLE算法流程（待实现）

1. 初始化场数据（速度、压力、温度）
2. 时间步循环：
   - 动量方程求解（u, v）
   - 压力修正方程求解
   - 速度和压力修正
   - 标量方程求解（温度等）
   - 检查收敛性
3. 输出结果

## 常用操作示例

```cpp
// 创建场
StructuredMesh mesh;
ScalarField temperature;
VectorField velocity;
BoundaryField bc(ncx, ncy, boundaryInfo, 4);
FluidPropertyField props;

// 访问数据
temperature(10, 20) = 300.0;
Velocity vel = velocity(10, 20);
FluidValues p = props(10, 20);

// Eigen求解
Eigen::MatrixXd A(3, 3);
Eigen::VectorXd b(3), x;
x = A.colPivHouseholderQr().solve(b);
```

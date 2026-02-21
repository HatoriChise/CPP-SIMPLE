cmake_minimum_required(VERSION 3.14)
project(CFD_Solve LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# 头文件目录（保留你的原有配置）
include_directories(${CMAKE_SOURCE_DIR})

# 关键修正1：显式设置CMP0167策略，消除Boost查找警告（适配新版CMake）
cmake_policy(SET CMP0167 NEW)

# 关键修正2：vcpkg库均加CONFIG关键字，指定Boost依赖组件（multi_array/math为项目必需）
find_package(fmt CONFIG REQUIRED)          # 显式CONFIG，适配vcpkg
find_package(Eigen3 CONFIG REQUIRED)       # 显式CONFIG，适配vcpkg
find_package(Boost CONFIG REQUIRED multi_array math)  # CONFIG+指定必需组件

# 项目源文件（保留你的原有配置）
set (SOURCES
    src/main.cpp
    src/test.cpp
    src/grid/structuredMesh.cpp
    src/field/scalarField.cpp
    src/field/vectorField.cpp
    src/field/boundaryField.cpp
    src/field/fluidProperties.cpp
)

# 生成可执行文件（保留你的原有配置）
add_executable(test_deps ${SOURCES})

# 关键修正3：添加Boost库链接（Boost::multi_array/Boost::math为组件对应的CMake目标）
# 保留原有fmt/Eigen3链接，按顺序追加Boost链接即可
target_link_libraries(test_deps PRIVATE
    fmt::fmt
    Eigen3::Eigen
    Boost::multi_array  # Boost多数组组件目标
    Boost::math         # Boost数学组件目标
)
target_compile_definitions(test_deps PRIVATE FMT_HEADER_ONLY=1)

# 编译选项（保留你的原有配置，无需修改）
if(CMAKE_BUILD_TYPE STREQUAL "Debug")
    target_compile_options(test_deps PRIVATE -O0 -g)
else()
    target_compile_options(test_deps PRIVATE -O2 -DNDEBUG)
endif()
cmake_minimum_required(VERSION 3.14)
project(CFD_Solve LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# headers
include_directories(${CMAKE_SOURCE_DIR})

find_package(fmt REQUIRED)
find_package(Eigen3 REQUIRED)
find_package(Boost REQUIRED)

set (SOURCES
    src/main.cpp
    src/test.cpp
    src/grid/structuredMesh.cpp
    src/field/scalarField.cpp
    src/field/vectorField.cpp
    src/field/boundaryField.cpp
    src/field/fluidProperties.cpp
)

add_executable(test_deps ${SOURCES})

target_link_libraries(test_deps PRIVATE fmt::fmt Eigen3::Eigen)

if(CMAKE_BUILD_TYPE STREQUAL "Debug")
    target_compile_options(test_deps PRIVATE -O0 -g)
else()
    target_compile_options(test_deps PRIVATE -O2 -DNDEBUG)
endif()
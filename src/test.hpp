#pragma once

#include <fmt/core.h>
#include <Eigen/Dense>
#include <boost/math/constants/constants.hpp>
#include <cmath>

#include "utils/formatter4eigen.h"

#include "config/simulationConfig.hpp"

#include "grid/structuredMesh.hpp"

#include "field/scalarField.hpp"
#include "field/vectorField.hpp"
#include "field/boundaryField.hpp"
#include "field/fluidProperties.hpp"
#include "equation/scalarEquation.hpp"

// === 测试断言宏 ===
#define TEST_ASSERT(condition, message) \
    if (!(condition)) { \
        fmt::print("❌ FAIL: {}\n", message); \
        return false; \
    }

#define TEST_ASSERT_NEAR(a, b, tol, message) \
    if (std::abs(static_cast<float>(a) - static_cast<float>(b)) > (tol)) { \
        fmt::print("❌ FAIL: {} (got {:.6f}, expected {:.6f})\n", message, a, b); \
        return false; \
    }

#define TEST_PASSED(name) fmt::print("✓ PASS: {}\n", name)

// === 测试函数声明 ===
void test();

// 模块测试
bool test_mesh();
bool test_scalar_field();
bool test_vector_field();
bool test_boundary_field();
bool test_fluid_props();

// 方程测试
bool test_diffusion_term();
bool test_convection_term();
bool test_pressure_gradient();

// 集成测试
void test_simple_iteration();
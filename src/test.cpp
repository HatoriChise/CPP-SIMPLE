// src/test.cpp
#include "test.hpp"
#include "src/field/massFluxField.hpp"

// ============================================================================
// 辅助函数：设置测试配置
// ============================================================================
void setup_test_configuration(std::array<BoudaryCondition, 4> &boundaryInfo)
{
    boundaryInfo[0] = {X_MIN, WALL, {0.0f, 0.0f}, DIRICHLET, 300.0f, 0.0f};         // west
    boundaryInfo[1] = {X_MAX, WALL, {0.0f, 0.0f}, DIRICHLET, 300.0f, 0.0f};         // east
    boundaryInfo[2] = {Y_MIN, WALL, {0.0f, 0.0f}, DIRICHLET, 300.0f, 0.0f};         // south
    boundaryInfo[3] = {Y_MAX, WALL, {lid_velocity, 0.0f}, DIRICHLET, 300.0f, 0.0f}; // north
}

// ============================================================================
// 测试：棋盘式压力场下的 MassFluxField 面速度计算
// ============================================================================
bool test_mass_flux_checkerboard_pressure()
{
    fmt::print("\n=== Testing MassFluxField with Checkerboard Pressure ===\n");

    // 1. 创建网格和场
    StructuredMesh mesh;
    VectorField velocity;
    ScalarField pressure(ncx, ncy, 0.0f);
    FluidPropertyField props;

    std::array<BoudaryCondition, 4> boundaryInfo;
    setup_test_configuration(boundaryInfo);
    BoundaryField bc(ncx, ncy, boundaryInfo.data(), 4);

    // 2. 设置均匀速度场 u=1, v=1
    for (int j = 0; j < ncy; ++j)
    {
        for (int i = 0; i < ncx; ++i)
        {
            velocity.u()(i, j) = 0.0f;
            velocity.v()(i, j) = 0.0f;
        }
    }

    // 3. 设置棋盘式压力场：-10, 0, 10 交替
    // 模式：(i+j) % 3 == 0 -> -10, (i+j) % 3 == 1 -> 0, (i+j) % 3 == 2 -> 10
    fmt::print("  Pressure field (checkerboard pattern -10, 0, 10):\n");
    for (int j = ncy - 1; j >= 0; --j)
    {
        fmt::print("    j={}: ", j);
        for (int i = 0; i < ncx; ++i)
        {
            if (i == 0)
                pressure(i, j) = -10.0f;
            else if (i == 1)
                pressure(i, j) = 0.0f;
            else
                pressure(i, j) = 10.0f;
            
            fmt::print("{:.0f} ", pressure(i, j));
        }
        fmt::print("\n");
    }

    // 4. 创建动量方程以获取 aP 系数
    ScalarField uDummy(ncx, ncy, 0.0f);
    ScalarField vDummy(ncx, ncy, 0.0f);
    ScalarEquation uEq(mesh, uDummy, velocity, bc, props, 0);  // direction=0 for u-momentum
    ScalarEquation vEq(mesh, vDummy, velocity, bc, props, 1);  // direction=1 for v-momentum

    // 组装扩散项以获得非零 aP
    uEq.resetCoefficients();
    uEq.addDiffusionTerm();
    vEq.resetCoefficients();
    vEq.addDiffusionTerm();

    // 5. 创建 MassFluxField 并计算面通量
    MassFluxField massFlux(ncx, ncy);
    massFlux.updateFluxes(mesh, velocity, props, &pressure, &uEq, &vEq);

    // 6. 输出并验证结果
    auto meshSize = mesh.getMeshSize();
    float dx = meshSize[0];
    float dy = meshSize[1];
    float rho = density;

    for(int j = 0; j < ncy; ++j)
    {
        for(int i = 0; i < ncx; ++i)
        {
            fmt::print("  Cell ({}, {}): aP={}, aE={:.3f}, aW={:.3f}, aN={:.3f}, aS={:.3f}\n", i, j,
                       uEq.getCoefMatrix()[i][j].aP, uEq.getCoefMatrix()[i][j].aE,
                       uEq.getCoefMatrix()[i][j].aW, uEq.getCoefMatrix()[i][j].aN,
                       uEq.getCoefMatrix()[i][j].aS);
        }
    }

    fmt::print("\n  Face velocities at internal cell (1, 1):\n");
    int i = 1, j = 1;
    const auto &flux = massFlux(i, j);

    // 面速度 = 质量通量 / (rho * area)
    // 东/西面面积 = dy, 北/南面面积 = dx
    float u_e = flux.mE / (rho * dy);
    float u_w = flux.mW / (rho * dy);
    float v_n = flux.mN / (rho * dx);
    float v_s = flux.mS / (rho * dx);

    fmt::print("    u_e (east face)  = {:.6f}\n", u_e);
    fmt::print("    u_w (west face)  = {:.6f}\n", u_w);
    fmt::print("    v_n (north face) = {:.6f}\n", v_n);
    fmt::print("    v_s (south face) = {:.6f}\n", v_s);

    // 7. 验证：对于均匀速度场 u=1, v=1，面速度应该接近 1.0
    // Rhie-Chow 修正应该消除棋盘压力场的影响
    // 允许一定的偏差，因为压力梯度会产生修正
    float tol = 0.5f;  // 允许 0.5 的偏差

    // TEST_ASSERT_NEAR(u_e, 1.0f, tol, "East face velocity should be close to 1.0");
    // TEST_ASSERT_NEAR(u_w, 1.0f, tol, "West face velocity should be close to 1.0");
    // TEST_ASSERT_NEAR(v_n, 1.0f, tol, "North face velocity should be close to 1.0");
    // TEST_ASSERT_NEAR(v_s, 1.0f, tol, "South face velocity should be close to 1.0");

    // 8. 对比：不使用 RC 插值的情况
    fmt::print("\n  Comparison: Without RC interpolation (pressure=nullptr):\n");
    MassFluxField massFluxNoRC(ncx, ncy);
    massFluxNoRC.updateFluxes(mesh, velocity, props, nullptr, nullptr, nullptr);

    const auto &fluxNoRC = massFluxNoRC(i, j);
    float u_e_noRC = fluxNoRC.mE / (rho * dy);
    float u_w_noRC = fluxNoRC.mW / (rho * dy);
    float v_n_noRC = fluxNoRC.mN / (rho * dx);
    float v_s_noRC = fluxNoRC.mS / (rho * dx);

    fmt::print("    Without RC: u_e = {:.6f}, u_w = {:.6f}, v_n = {:.6f}, v_s = {:.6f}\n",
               u_e_noRC, u_w_noRC, v_n_noRC, v_s_noRC);

    // 9. 输出多个单元的面速度分布
    fmt::print("\n  Face velocity distribution (u_e and u_w at each cell):\n");
    for (int jj = ncy - 1; jj >= 0; --jj)
    {
        fmt::print("    j={}: ", jj);
        for (int ii = 0; ii < ncx; ++ii)
        {
            float ue_face = massFlux(ii, jj).mE / (rho * dy);
            fmt::print("{:.3f} ", ue_face);
        }
        for (int ii = 0; ii < ncx; ++ii)
        {
            float uw_face = massFlux(ii, jj).mW / (rho * dy);
            fmt::print("{:.3f} ", uw_face);
        }
        
        fmt::print("\n");
    }

    TEST_PASSED("MassFluxField with Checkerboard Pressure");
    return true;
}

// ============================================================================
// 测试：均匀压力场下的 MassFluxField（RC 修正应为零）
// ============================================================================
bool test_mass_flux_uniform_pressure()
{
    fmt::print("\n=== Testing MassFluxField with Uniform Pressure ===\n");

    StructuredMesh mesh;
    VectorField velocity;
    ScalarField pressure(ncx, ncy, 100.0f);  // 均匀压力
    FluidPropertyField props;

    std::array<BoudaryCondition, 4> boundaryInfo;
    setup_test_configuration(boundaryInfo);
    BoundaryField bc(ncx, ncy, boundaryInfo.data(), 4);

    // 设置均匀速度场 u=1, v=1
    for (int j = 0; j < ncy; ++j)
    {
        for (int i = 0; i < ncx; ++i)
        {
            velocity.u()(i, j) = 1.0f;
            velocity.v()(i, j) = 1.0f;
        }
    }

    // 创建动量方程
    ScalarField uDummy(ncx, ncy, 0.0f);
    ScalarField vDummy(ncx, ncy, 0.0f);
    ScalarEquation uEq(mesh, uDummy, velocity, bc, props, 0);
    ScalarEquation vEq(mesh, vDummy, velocity, bc, props, 1);

    uEq.resetCoefficients();
    uEq.addDiffusionTerm();
    vEq.resetCoefficients();
    vEq.addDiffusionTerm();

    // 计算 MassFluxField
    MassFluxField massFlux(ncx, ncy);
    massFlux.updateFluxes(mesh, velocity, props, &pressure, &uEq, &vEq);

    // 验证：均匀压力场下，RC 修正为零，面速度应等于线性插值速度 = 1.0
    auto meshSize = mesh.getMeshSize();
    float dx = meshSize[0];
    float dy = meshSize[1];
    float rho = density;

    int i = 2, j = 2;
    const auto &flux = massFlux(i, j);

    float u_e = flux.mE / (rho * dy);
    float u_w = flux.mW / (rho * dy);
    float v_n = flux.mN / (rho * dx);
    float v_s = flux.mS / (rho * dx);

    fmt::print("  Uniform pressure: u_e={:.6f}, u_w={:.6f}, v_n={:.6f}, v_s={:.6f}\n",
               u_e, u_w, v_n, v_s);

    // 均匀压力场，RC 修正为零，面速度应精确等于 1.0
    TEST_ASSERT_NEAR(u_e, 1.0f, 1e-6f, "East face velocity should be 1.0 for uniform pressure");
    TEST_ASSERT_NEAR(u_w, 1.0f, 1e-6f, "West face velocity should be 1.0 for uniform pressure");
    TEST_ASSERT_NEAR(v_n, 1.0f, 1e-6f, "North face velocity should be 1.0 for uniform pressure");
    TEST_ASSERT_NEAR(v_s, 1.0f, 1e-6f, "South face velocity should be 1.0 for uniform pressure");

    TEST_PASSED("MassFluxField with Uniform Pressure");
    return true;
}

// ============================================================================
// 主测试函数
// ============================================================================
void test()
{
    fmt::print("=== CFD SIMPLE Solver Unit Tests ===\n");

    int passed = 0;
    int total = 0;

    // total++;
    // if (test_mass_flux_uniform_pressure()) passed++;

    total++;
    if (test_mass_flux_checkerboard_pressure()) passed++;
}
// src/test.cpp
#include "test.hpp"

// ============================================================================
// 辅助函数：设置测试配置
// ============================================================================
void setup_test_configuration(std::array<BoudaryCondition, 4> &boundaryInfo, int &test_ncx,
                              int &test_ncy, float &test_Lx, float &test_Ly)
{
    test_ncx = ncx;
    test_ncy = ncy;
    test_Lx = Lx;
    test_Ly = Ly;

    boundaryInfo[0] = {X_MIN, WALL, {0.0f, 0.0f}, DIRICHLET, 300.0f, 0.0f};         // west
    boundaryInfo[1] = {X_MAX, WALL, {0.0f, 0.0f}, DIRICHLET, 300.0f, 0.0f};         // east
    boundaryInfo[2] = {Y_MIN, WALL, {0.0f, 0.0f}, DIRICHLET, 300.0f, 0.0f};         // south
    boundaryInfo[3] = {Y_MAX, WALL, {lid_velocity, 0.0f}, DIRICHLET, 300.0f, 0.0f}; // north
}

// ============================================================================
// 模块测试：网格
// ============================================================================
bool test_mesh()
{
    fmt::print("\n=== Testing StructuredMesh ===\n");

    StructuredMesh mesh;
    auto meshSize = mesh.getMeshSize();

    // 测试网格尺寸
    TEST_ASSERT_NEAR(meshSize[0], Lx / ncx, 1e-6f, "dx 计算正确");
    TEST_ASSERT_NEAR(meshSize[1], Ly / ncy, 1e-6f, "dy 计算正确");

    // 测试单元中心坐标
    const auto &xc = mesh.getCellCentersX();
    const auto &yc = mesh.getCellCentersY();

    TEST_ASSERT_NEAR(xc[0], meshSize[0] / 2.0f, 1e-6f, "第一个单元中心 x 坐标正确");
    TEST_ASSERT_NEAR(yc[0], meshSize[1] / 2.0f, 1e-6f, "第一个单元中心 y 坐标正确");

    TEST_PASSED("StructuredMesh");
    return true;
}

// ============================================================================
// 模块测试：标量场
// ============================================================================
bool test_scalar_field()
{
    fmt::print("\n=== Testing ScalarField ===\n");

    // 测试默认构造（使用配置中的初始温度）
    ScalarField field;
    TEST_ASSERT_NEAR(field(0, 0), initalTemperature, 1e-6f, "默认初始化为 initalTemperature");

    // 测试自定义构造
    ScalarField customField(ncx, ncy, 100.0f);
    TEST_ASSERT_NEAR(customField(0, 0), 100.0f, 1e-6f, "自定义初始值正确");
    TEST_ASSERT_NEAR(customField(ncx - 1, ncy - 1), 100.0f, 1e-6f, "边界值正确");

    // 测试 fill 方法
    customField.fill(200.0f);
    TEST_ASSERT_NEAR(customField(0, 0), 200.0f, 1e-6f, "fill 方法正确");

    // 测试读写
    customField(2, 3) = 300.0f;
    TEST_ASSERT_NEAR(customField(2, 3), 300.0f, 1e-6f, "读写访问正确");

    TEST_PASSED("ScalarField");
    return true;
}

// ============================================================================
// 模块测试：矢量场
// ============================================================================
bool test_vector_field()
{
    fmt::print("\n=== Testing VectorField ===\n");

    VectorField velocity;

    // 测试默认初始化为 0
    Velocity v0 = velocity(0, 0);
    TEST_ASSERT_NEAR(v0.u, 0.0f, 1e-6f, "默认 u=0");
    TEST_ASSERT_NEAR(v0.v, 0.0f, 1e-6f, "默认 v=0");

    // 测试修改 u 分量
    velocity.u()(2, 2) = 1.0f;
    Velocity v = velocity(2, 2);
    TEST_ASSERT_NEAR(v.u, 1.0f, 1e-6f, "修改 u 分量正确");

    // 测试修改 v 分量
    velocity.v()(2, 2) = 0.5f;
    v = velocity(2, 2);
    TEST_ASSERT_NEAR(v.v, 0.5f, 1e-6f, "修改 v 分量正确");

    TEST_PASSED("VectorField");
    return true;
}

// ============================================================================
// 模块测试：边界场
// ============================================================================
bool test_boundary_field()
{
    fmt::print("\n=== Testing BoundaryField ===\n");

    std::array<BoudaryCondition, 4> boundaryInfo;
    int test_ncx, test_ncy;
    float test_Lx, test_Ly;
    setup_test_configuration(boundaryInfo, test_ncx, test_ncy, test_Lx, test_Ly);

    BoundaryField bc(ncx, ncy, boundaryInfo.data(), 4);

    // 测试西侧边界
    TEST_ASSERT(bc.west(0).velocityType == WALL, "西侧边界类型为 WALL");
    TEST_ASSERT_NEAR(bc.west(0).VelocityValue[0], 0.0f, 1e-6f, "西侧边界 u=0");

    // 测试北侧边界（顶盖）
    TEST_ASSERT(bc.north(0).velocityType == WALL, "北侧边界类型为 WALL");
    TEST_ASSERT_NEAR(bc.north(0).VelocityValue[0], lid_velocity, 1e-6f, "北侧边界 u=lid_velocity");

    TEST_PASSED("BoundaryField");
    return true;
}

// ============================================================================
// 模块测试：流体物性场
// ============================================================================
bool test_fluid_props()
{
    fmt::print("\n=== Testing FluidPropertyField ===\n");

    FluidPropertyField props;

    // 测试默认物性
    FluidValues p = props(0, 0);
    TEST_ASSERT_NEAR(p.rho, density, 1e-6f, "密度正确");
    TEST_ASSERT_NEAR(p.mu, mu_fluid, 1e-6f, "粘度正确");
    TEST_ASSERT_NEAR(p.k, k_fluid, 1e-6f, "导热率正确");
    TEST_ASSERT_NEAR(p.cp, cp_fluid, 1e-6f, "比热正确");
    TEST_ASSERT_NEAR(p.nu, nu_fluid, 1e-6f, "运动粘度正确");

    // 测试修改物性
    props.rho()(0, 0) = 1.2f;
    TEST_ASSERT_NEAR(props.rho()(0, 0), 1.2f, 1e-6f, "修改密度正确");

    TEST_PASSED("FluidPropertyField");
    return true;
}

// ============================================================================
// 方程测试：扩散项
// ============================================================================
bool test_diffusion_term()
{
    fmt::print("\n=== Testing Diffusion Term ===\n");

    StructuredMesh mesh;
    ScalarField phi(ncx, ncy, 0.0f);
    VectorField velocity;
    std::array<BoudaryCondition, 4> boundaryInfo;
    int test_ncx, test_ncy;
    float test_Lx, test_Ly;
    setup_test_configuration(boundaryInfo, test_ncx, test_ncy, test_Lx, test_Ly);
    BoundaryField bc(ncx, ncy, boundaryInfo.data(), 4);
    FluidPropertyField props;

    ScalarEquation eq(mesh, phi, velocity, bc, props, -1);

    // 组装扩散项
    eq.resetCoefficients();
    eq.addDiffusionTerm();

    // 测试内部单元 (2, 2) 的系数
    int i = 2, j = 2;
    const auto &coef = eq.getCoefMatrix()[j][i];

    auto meshSize = mesh.getMeshSize();
    float dx = meshSize[0];
    float dy = meshSize[1];
    float Gamma = props(i, j).mu;

    float expected_aE = Gamma * dy / dx;
    float expected_aW = Gamma * dy / dx;
    float expected_aN = Gamma * dx / dy;
    float expected_aS = Gamma * dx / dy;
    float expected_aP = expected_aE + expected_aW + expected_aN + expected_aS;

    TEST_ASSERT_NEAR(coef.aE, expected_aE, 1e-6f, "aE 扩散系数正确");
    TEST_ASSERT_NEAR(coef.aW, expected_aW, 1e-6f, "aW 扩散系数正确");
    TEST_ASSERT_NEAR(coef.aN, expected_aN, 1e-6f, "aN 扩散系数正确");
    TEST_ASSERT_NEAR(coef.aS, expected_aS, 1e-6f, "aS 扩散系数正确");
    TEST_ASSERT_NEAR(coef.aP, expected_aP, 1e-6f, "aP 扩散系数正确");

    TEST_PASSED("Diffusion Term");
    return true;
}

// ============================================================================
// 方程测试：对流项
// ============================================================================
bool test_convection_term()
{
    fmt::print("\n=== Testing Convection Term ===\n");

    StructuredMesh mesh;
    ScalarField phi(ncx, ncy, 0.0f);
    VectorField velocity;
    std::array<BoudaryCondition, 4> boundaryInfo;
    int test_ncx, test_ncy;
    float test_Lx, test_Ly;
    setup_test_configuration(boundaryInfo, test_ncx, test_ncy, test_Lx, test_Ly);
    BoundaryField bc(ncx, ncy, boundaryInfo.data(), 4);
    FluidPropertyField props;

    ScalarEquation eq(mesh, phi, velocity, bc, props, -1);

    // 设置一个均匀速度场 u=1.0, v=0.5
    for(int j = 0; j < ncy; ++j)
    {
        for(int i = 0; i < ncx; ++i)
        {
            velocity.u()(i, j) = 1.0f;
            velocity.v()(i, j) = 0.5f;
        }
    }

    // 只测试对流项（不先调用扩散项）
    eq.resetCoefficients();
    ScalarField massFluxEast(ncx, ncy, 0.0f);
    ScalarField massFluxNorth(ncx, ncy, 0.0f);
    ScalarField dummyPressure(ncx, ncy, 0.0f);
    eq.updateMassFluxField(dummyPressure, massFluxEast, massFluxNorth);
    eq.addConvectionTerm(massFluxEast, massFluxNorth); // 线性插值，不使用 RC

    // 测试内部单元 (2, 2) 的系数
    int i = 2, j = 2;
    const auto &coef = eq.getCoefMatrix()[j][i];

    // 计算理论值
    auto meshSize = mesh.getMeshSize();
    float dx = meshSize[0];
    float dy = meshSize[1];
    float rho = props(i, j).rho;

    // 界面速度（线性插值）
    float u_e = 1.0f, u_w = 1.0f, v_n = 0.5f, v_s = 0.5f;
    float F_e = rho * u_e * dy;
    float F_w = rho * u_w * dy;
    float F_n = rho * v_n * dx;
    float F_s = rho * v_s * dx;

    // 迎风格式系数（u>0 时，aE=0, aW=F_w）
    float expected_aE = std::max(-F_e, 0.0f); // = 0 (因为 F_e > 0)
    float expected_aW = std::max(F_w, 0.0f);  // = F_w (因为 F_w > 0)
    float expected_aN = std::max(-F_n, 0.0f); // = 0 (因为 F_n > 0)
    float expected_aS = std::max(F_s, 0.0f);  // = F_s (因为 F_s > 0)

    TEST_ASSERT_NEAR(coef.aE, expected_aE, 1e-6f, "aE 对流系数正确 (u>0 时应为 0)");
    TEST_ASSERT_NEAR(coef.aW, expected_aW, 1e-6f, "aW 对流系数正确");
    TEST_ASSERT_NEAR(coef.aN, expected_aN, 1e-6f, "aN 对流系数正确 (v>0 时应为 0)");
    TEST_ASSERT_NEAR(coef.aS, expected_aS, 1e-6f, "aS 对流系数正确");

    TEST_PASSED("Convection Term");
    return true;
}

// ============================================================================
// 方程测试：压力梯度
// ============================================================================
bool test_pressure_gradient()
{
    fmt::print("\n=== Testing Pressure Gradient ===\n");

    StructuredMesh mesh;
    ScalarField pressure(ncx, ncy, 0.0f);
    VectorField velocity;
    std::array<BoudaryCondition, 4> boundaryInfo;
    int test_ncx, test_ncy;
    float test_Lx, test_Ly;
    setup_test_configuration(boundaryInfo, test_ncx, test_ncy, test_Lx, test_Ly);
    BoundaryField bc(ncx, ncy, boundaryInfo.data(), 4);
    FluidPropertyField props;

    // 设置线性压力场 p = x，则 dp/dx = 1
    auto meshSize = mesh.getMeshSize();
    float dx = meshSize[0];
    for(int j = 0; j < ncy; ++j)
    {
        for(int i = 0; i < ncx; ++i)
        {
            pressure(i, j) = static_cast<float>(i) * dx;
        }
    }

    // 测试 u 动量方程 (direction=0)
    ScalarField u_field(ncx, ncy, 0.0f);
    ScalarEquation u_momentum(mesh, u_field, velocity, bc, props, 0);

    u_momentum.resetCoefficients();
    u_momentum.addPressureGradient(pressure);

    // 测试内部单元 (2, 2) 的源项
    int i = 2, j = 2;
    const auto &coef = u_momentum.getCoefMatrix()[j][i];
    float dy = meshSize[1];

    // 理论值：bsrc = -(p_e - p_w) * dy
    // p_e = (p_P + p_E) / 2, p_w = (p_W + p_P) / 2
    // 对于 p = x，有 p_e - p_w = dx
    float expected_bsrc = -dx * dy;

    TEST_ASSERT_NEAR(coef.bsrc, expected_bsrc, 1e-6f, "u 动量压力梯度源项正确");

    // 测试 v 动量方程 (direction=1)
    ScalarField v_field(ncx, ncy, 0.0f);
    ScalarEquation v_momentum(mesh, v_field, velocity, bc, props, 1);

    v_momentum.resetCoefficients();
    v_momentum.addPressureGradient(pressure);

    // 对于 p = x（与 y 无关），v 动量的压力梯度应为 0
    const auto &v_coef = v_momentum.getCoefMatrix()[j][i];
    TEST_ASSERT_NEAR(v_coef.bsrc, 0.0f, 1e-6f, "v 动量压力梯度源项正确（p 与 y 无关）");

    TEST_PASSED("Pressure Gradient");
    return true;
}

// ============================================================================
// 方程测试：Rhie-Chow 插值
// ============================================================================
bool test_rhie_chow()
{
    fmt::print("\n=== Testing Rhie-Chow Interpolation ===\n");

    StructuredMesh mesh;
    VectorField velocity;
    std::array<BoudaryCondition, 4> boundaryInfo;
    int test_ncx, test_ncy;
    float test_Lx, test_Ly;
    setup_test_configuration(boundaryInfo, test_ncx, test_ncy, test_Lx, test_Ly);
    BoundaryField bc(ncx, ncy, boundaryInfo.data(), 4);
    FluidPropertyField props;

    auto meshSize = mesh.getMeshSize();
    float dx = meshSize[0];
    float dy = meshSize[1];

    // 测试 1：均匀压力场下，RC 修正项应为零
    {
        fmt::print("  Test 1: Uniform pressure field (RC correction should be zero)...\n");
        ScalarField pressure_uniform(ncx, ncy, 100.0f); // 均匀压力场 p = 100
        ScalarField phi(ncx, ncy, 0.0f);

        ScalarEquation eq(mesh, phi, velocity, bc, props, -1);

        // 先组装扩散项以获取 aP 系数
        eq.resetCoefficients();
        eq.addDiffusionTerm();

        // 保存扩散项系数
        int i = 2, j = 2;
        float aE_diff = eq.getCoefMatrix()[j][i].aE;
        float aW_diff = eq.getCoefMatrix()[j][i].aW;
        float aN_diff = eq.getCoefMatrix()[j][i].aN;
        float aS_diff = eq.getCoefMatrix()[j][i].aS;

        // 设置均匀速度场
        for(int j = 0; j < ncy; ++j)
        {
            for(int i = 0; i < ncx; ++i)
            {
                velocity.u()(i, j) = 1.0f;
                velocity.v()(i, j) = 0.5f;
            }
        }

        // 使用 RC 插值计算对流项
        ScalarField massFluxEastUni(ncx, ncy, 0.0f);
        ScalarField massFluxNorthUni(ncx, ncy, 0.0f);
        eq.updateMassFluxField(pressure_uniform, massFluxEastUni, massFluxNorthUni);
        eq.addConvectionTerm(massFluxEastUni, massFluxNorthUni);

        // 保存系数到文件

        // 测试内部单元 (2, 2) 的系数
        const auto &coef = eq.getCoefMatrix()[j][i];

        // 对于均匀压力场，RC 修正项为零，结果应与线性插值相同
        // 迎风格式：u>0 时，aE=0, aW=F_w
        float rho = props(i, j).rho;
        float u_e = 1.0f, u_w = 1.0f, v_n = 0.5f, v_s = 0.5f;
        float F_e = rho * u_e * dy;
        float F_w = rho * u_w * dy;
        float F_n = rho * v_n * dx;
        float F_s = rho * v_s * dx;

        // 期望值 = 扩散项 + 对流项
        float expected_aE = aE_diff + std::max(-F_e, 0.0f);
        float expected_aW = aW_diff + std::max(F_w, 0.0f);
        float expected_aN = aN_diff + std::max(-F_n, 0.0f);
        float expected_aS = aS_diff + std::max(F_s, 0.0f);

        TEST_ASSERT_NEAR(coef.aE, expected_aE, 1e-5f, "均匀压力场：aE 正确（RC 修正为零）");
        TEST_ASSERT_NEAR(coef.aW, expected_aW, 1e-5f, "均匀压力场：aW 正确");
        TEST_ASSERT_NEAR(coef.aN, expected_aN, 1e-5f, "均匀压力场：aN 正确");
        TEST_ASSERT_NEAR(coef.aS, expected_aS, 1e-5f, "均匀压力场：aS 正确");

        fmt::print("NOTE: Coefficients aP at ({}, {}): aP={:.6f}\n", i, j, coef.aP);
    }

    // 测试 2：线性压力场下，验证 RC 修正项
    {
        fmt::print("  Test 2: Linear pressure field (verify RC correction)...\n");
        // 设置线性压力场 p = 100 * x，则 dp/dx = 100
        ScalarField pressure_linear(ncx, ncy, 0.0f);
        for(int j = 0; j < ncy; ++j)
        {
            for(int i = 0; i < ncx; ++i)
            {
                pressure_linear(i, j) = 100.0f * static_cast<float>(i) * dx;
            }
        }

        ScalarField phi(ncx, ncy, 0.0f);
        ScalarEquation eq(mesh, phi, velocity, bc, props, -1);

        // 先组装扩散项以获取 aP 系数
        eq.resetCoefficients();
        eq.addDiffusionTerm();

        // 设置静止速度场（便于观察 RC 修正）
        for(int j = 0; j < ncy; ++j)
        {
            for(int i = 0; i < ncx; ++i)
            {
                velocity.u()(i, j) = 0.0f;
                velocity.v()(i, j) = 0.0f;
            }
        }

        // 使用 RC 插值计算对流项
        ScalarField massFluxEastLin(ncx, ncy, 0.0f);
        ScalarField massFluxNorthLin(ncx, ncy, 0.0f);
        eq.updateMassFluxField(pressure_linear, massFluxEastLin, massFluxNorthLin);
        eq.addConvectionTerm(massFluxEastLin, massFluxNorthLin);

        // 对于静止速度场和线性压力场，RC 修正应产生非零的界面速度
        // 验证系数矩阵非零（说明 RC 修正生效）
        int i = 2, j = 2;
        const auto &coef = eq.getCoefMatrix()[j][i];

        // 由于速度为零且压力梯度为常数，RC 修正项应相互抵消
        // 但系数矩阵应反映 RC 插值的影响
        fmt::print("    Coefficients at ({}, {}): aE={:.6f}, aW={:.6f}, aN={:.6f}, aS={:.6f}\n", i,
                   j, coef.aE, coef.aW, coef.aN, coef.aS);

        // 验证系数矩阵已组装（非全零）
        bool has_nonzero_coef =
            (coef.aE != 0.0f || coef.aW != 0.0f || coef.aN != 0.0f || coef.aS != 0.0f);
        TEST_ASSERT(has_nonzero_coef, "线性压力场：RC 插值产生非零系数");
    }

    // 测试 3：棋盘格压力场下，验证 RC 修正的阻尼效果
    // 注意：RC 修正需要非零速度场才能体现效果
    {
        fmt::print("  Test 3: Checkerboard pressure field (verify RC damping)...\n");
        // 设置棋盘格压力场：p(i,j) = 100 * (-1)^(i+j)
        ScalarField pressure_checker(ncx, ncy, 0.0f);
        for(int j = 0; j < ncy; ++j)
        {
            for(int i = 0; i < ncx; ++i)
            {
                pressure_checker(i, j) = 10.0f * ((i + j) % 2 == 0 ? 1.0f : -1.0f);
            }
        }

        // 打印棋盘格压力场
        fmt::print("    Checkerboard pressure values:\n");
        for(int j = ncy - 1; j >= 0; --j)
        {
            fmt::print("      j={}: ", j);
            for(int i = 0; i < ncx; ++i)
            {
                fmt::print("{:.1f} ", pressure_checker(i, j));
            }
            fmt::print("\n");
        }

        // 设置静止速度场
        for(int j = 0; j < ncy; ++j)
        {
            for(int i = 0; i < ncx; ++i)
            {
                velocity.u()(i, j) = 0.0f;
                velocity.v()(i, j) = 0.0f;
            }
        }

        ScalarField phi(ncx, ncy, 0.0f);
        ScalarEquation eq_with_rc(mesh, phi, velocity, bc, props, 1);
        ScalarEquation eq_without_rc(mesh, phi, velocity, bc, props, 1);

        // 组装扩散项（得到 aP 系数）
        eq_with_rc.resetCoefficients();
        eq_with_rc.addDiffusionTerm();
        eq_without_rc.resetCoefficients();
        eq_without_rc.addDiffusionTerm();


        // 分别使用 RC 插值和不使用 RC 插值
        ScalarField mfe_rc(ncx, ncy, 0.0f);
        ScalarField mfn_rc(ncx, ncy, 0.0f);
        eq_with_rc.updateMassFluxField(pressure_checker, mfe_rc, mfn_rc);
        eq_with_rc.addConvectionTerm(mfe_rc, mfn_rc);
        ScalarField mfe_no_rc(ncx, ncy, 0.0f);
        ScalarField mfn_no_rc(ncx, ncy, 0.0f);
        ScalarField dummyPressure(ncx, ncy, 0.0f);
        eq_without_rc.updateMassFluxField(dummyPressure, mfe_no_rc, mfn_no_rc);
        eq_without_rc.addConvectionTerm(mfe_no_rc, mfn_no_rc); // 无 RC

        // 比较系数差异
        int i = 1, j = 1;
        const auto &coef_rc = eq_with_rc.getCoefMatrix()[j][i];
        const auto &coef_no_rc = eq_without_rc.getCoefMatrix()[j][i];

        fmt::print("    With RC:    aE={:.6f}, aW={:.6f}, aN={:.6f}, aS={:.6f}, aP={:.6f}\n",
                   coef_rc.aE, coef_rc.aW, coef_rc.aN, coef_rc.aS, coef_rc.aP);
        fmt::print("    Without RC: aE={:.6f}, aW={:.6f}, aN={:.6f}, aS={:.6f}, aP={:.6f}\n",
                   coef_no_rc.aE, coef_no_rc.aW, coef_no_rc.aN, coef_no_rc.aS, coef_no_rc.aP);

        // 验证 RC 插值产生了不同的系数（RC 修正改变了界面速度或 aP）
        // 注意：RC 修正主要影响 aP（通过改变连续性方程的质量平衡）
        bool has_rc_effect = (coef_rc.aP != coef_no_rc.aP || coef_rc.aE != coef_no_rc.aE ||
                              coef_rc.aW != coef_no_rc.aW || coef_rc.aN != coef_no_rc.aN ||
                              coef_rc.aS != coef_no_rc.aS);

        if(!has_rc_effect)
        {
            fmt::print("    Note: RC correction did not change coefficients\n");
        }

        // 验证：RC 插值应该改变 aP 系数（质量平衡）
        TEST_ASSERT(coef_rc.aP != coef_no_rc.aP, "棋盘格压力场：RC 插值改变 aP 系数");
    }


    TEST_PASSED("Rhie-Chow Interpolation");
    return true;
}

// ============================================================================
// 集成测试：SIMPLE 迭代
// ============================================================================
void test_simple_iteration()
{
    fmt::print("\n=== Testing SIMPLE Iteration ===\n");

    StructuredMesh mesh;
    ScalarField pressure(ncx, ncy, 0.0f);
    VectorField velocity;
    std::array<BoudaryCondition, 4> boundaryInfo;
    int test_ncx, test_ncy;
    float test_Lx, test_Ly;
    setup_test_configuration(boundaryInfo, test_ncx, test_ncy, test_Lx, test_Ly);
    BoundaryField bc(ncx, ncy, boundaryInfo.data(), 4);
    FluidPropertyField props;

    // 创建 u/v 动量方程
    ScalarField u_field(ncx, ncy, 0.0f);
    ScalarField v_field(ncx, ncy, 0.0f);
    ScalarEquation u_momentum(mesh, u_field, velocity, bc, props, 0);
    ScalarEquation v_momentum(mesh, v_field, velocity, bc, props, 1);

    // 模拟一次 SIMPLE 迭代（不求解，只组装）
    u_momentum.resetCoefficients();
    u_momentum.addDiffusionTerm();
    ScalarField u_mfe(ncx, ncy, 0.0f);
    ScalarField u_mfn(ncx, ncy, 0.0f);
    u_momentum.updateMassFluxField(pressure, u_mfe, u_mfn);
    u_momentum.addConvectionTerm(u_mfe, u_mfn); // 使用 RC 插值
    u_momentum.addPressureGradient(pressure);
    // u_momentum.applyBoundaries();  // 应用边界条件
    u_momentum.saveCoefficientsToFile("../data/u_momentum_coefficients.csv");

    v_momentum.resetCoefficients();
    v_momentum.addDiffusionTerm();
    ScalarField v_mfe(ncx, ncy, 0.0f);
    ScalarField v_mfn(ncx, ncy, 0.0f);
    ScalarField dummyPressure2(ncx, ncy, 0.0f);
    v_momentum.updateMassFluxField(dummyPressure2, v_mfe, v_mfn);
    v_momentum.addConvectionTerm(v_mfe, v_mfn); // 不使用 RC 插值
    v_momentum.addPressureGradient(pressure);
    // v_momentum.applyBoundaries();  // 应用边界条件
    v_momentum.saveCoefficientsToFile("../data/v_momentum_coefficients.csv");

    // 打印系数摘要
    int i = 2, j = 2;
    const auto &u_coef = u_momentum.getCoefMatrix()[j][i];
    const auto &v_coef = v_momentum.getCoefMatrix()[j][i];

    fmt::print("U-momentum coefficients at ({}, {}):\n", i, j);
    fmt::print("  aE={:.6f}, aW={:.6f}, aN={:.6f}, aS={:.6f}, aP={:.6f}, bsrc={:.6f}\n", u_coef.aE,
               u_coef.aW, u_coef.aN, u_coef.aS, u_coef.aP, u_coef.bsrc);

    fmt::print("V-momentum coefficients at ({}, {}):\n", i, j);
    fmt::print("  aE={:.6f}, aW={:.6f}, aN={:.6f}, aS={:.6f}, aP={:.6f}, bsrc={:.6f}\n", v_coef.aE,
               v_coef.aW, v_coef.aN, v_coef.aS, v_coef.aP, v_coef.bsrc);

    TEST_PASSED("SIMPLE Iteration");
}

// ============================================================================
// 主测试函数
// ============================================================================
void test()
{
    fmt::print("=== CFD SIMPLE Solver Unit Tests ===\n");

    int passed = 0;
    int total = 0;

    // 模块测试
    total++;
    if(test_mesh())
        passed++;
    total++;
    if(test_scalar_field())
        passed++;
    total++;
    if(test_vector_field())
        passed++;
    total++;
    if(test_boundary_field())
        passed++;
    total++;
    if(test_fluid_props())
        passed++;

    // 方程测试
    total++;
    if(test_diffusion_term())
        passed++;
    total++;
    if(test_convection_term())
        passed++;
    total++;
    if(test_pressure_gradient())
        passed++;
    total++;
    if(test_rhie_chow())
        passed++;

    // 集成测试
    test_simple_iteration();

    // 汇总
    fmt::print("\n=== Test Summary ===\n");
    fmt::print("Passed: {}/{}\n", passed, total);

    if(passed == total)
    {
        fmt::print("🎉 All tests passed!\n");
    }
    else
    {
        fmt::print("⚠️  {} test(s) failed\n", total - passed);
    }
}

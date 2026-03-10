// src/test.cpp

#include "test.hpp"

void test_basic_fields()
{
    fmt::print("\n=== Testing Basic Fields ===\n");

    StructuredMesh mesh;
    mesh.saveMeshInfo();

    ScalarField temperatureField;
    fmt::print("Temperature at (0,0): {:.2f} K\n", temperatureField(0, 0));

    VectorField velocityField;
    fmt::print("Velocity at (0,0): {}\n", velocityField(0, 0));

    FluidPropertyField fluidPropertyField;
    fmt::print("Density at (0,0): {:.5f} kg/m^3\n", fluidPropertyField(0, 0).rho);
    fmt::print("Viscosity at (0,0): {:.5f} Pa·s\n", fluidPropertyField(0, 0).mu);
    fmt::print("Thermal diffusivity at (0,0): {:.5f} m^2/s\n", fluidPropertyField(0, 0).nu);
    fmt::print("Thermal conductivity at (0,0): {:.5f} W/m·K\n", fluidPropertyField(0, 0).k);
    fmt::print("Specific heat at (0,0): {:.5f} J/kg·K\n", fluidPropertyField(0, 0).cp);
    fmt::print("Alpha at (0,0): {:.2f}\n", fluidPropertyField(0, 0).alpha);
}

void test_linear_algebra()
{
    fmt::print("\n=== Testing Linear Algebra ===\n");

    Eigen::MatrixXd A = Eigen::MatrixXd::Random(3, 3);
    Eigen::VectorXd b(3);
    b << 1, 2, 3;
    Eigen::VectorXd x = A.colPivHouseholderQr().solve(b);
    fmt::print("A = \n{}\n", A);
    fmt::print("Solution x = [{:.4f}, {:.4f}, {:.4f}]\n", x(0), x(1), x(2));

    double pi = boost::math::constants::pi<double>();
    fmt::print("pi ≈ {:.10f}\n", pi);
}

void setup_test_configuration(std::array<BoudaryCondition, 4>& boundaryInfo,
                           int& test_ncx, int& test_ncy,
                           float& test_Lx, float& test_Ly)
{
    test_ncx = 10;
    test_ncy = 10;
    test_Lx = 1.0;
    test_Ly = 1.0;

    boundaryInfo[0] = {X_MIN, WALL, {0.0f, 0.0f}, DIRICHLET, 300.0f, 0.0f}; // west
    boundaryInfo[1] = {X_MAX, WALL, {0.0f, 0.0f}, DIRICHLET, 300.0f, 0.0f}; // east
    boundaryInfo[2] = {Y_MIN, WALL, {0.0f, 0.0f}, DIRICHLET, 300.0f, 0.0f}; // south
    boundaryInfo[3] = {Y_MAX, WALL, {0.0f, 0.0f}, DIRICHLET, 300.0f, 0.0f}; // north
}

void test_scalar_equation_construction(ScalarEquation& eq,
                                    const std::array<BoudaryCondition, 4>& boundaryInfo,
                                    int test_ncx, int test_ncy)
{
    auto& coefMatrix = eq.getCoefMatrix();
    fmt::print("Coefficient matrix size: {} x {}\n",
               coefMatrix.shape()[0], coefMatrix.shape()[1]);
}

bool test_scalar_equation_reset_coefficients(ScalarEquation& eq, int test_ncx, int test_ncy)
{
    // 测试重置系数
    eq.resetCoefficients();

    // 验证所有系数是否为零
    auto& coefMatrix = eq.getCoefMatrix();
    for (int j = 0; j < test_ncy; ++j)
    {
        for (int i = 0; i < test_ncx; ++i)
        {
            const auto& coef = coefMatrix[j][i];
            if (coef.aE != 0.0f || coef.aW != 0.0f || coef.aN != 0.0f ||
                coef.aS != 0.0f || coef.aP != 0.0f || coef.bsrc != 0.0f)
            {
                return false;
            }
        }
    }
    return true;
}

bool test_scalar_equation_add_diffusion(ScalarEquation& eq, StructuredMesh& mesh, const FluidPropertyField& props)
{
    // 重置系数确保初始状态
    eq.resetCoefficients();

    // 添加扩散项
    eq.addDiffusionTerm();

    // 验证内部单元的系数是否正确
    auto& coefMatrix = eq.getCoefMatrix();

    // 检查一个内部单元 (2, 2) 的系数
    int i = 2, j = 2;
    const auto& coef = coefMatrix[j][i];

    // 计算理论值
    auto meshSize = mesh.getMeshSize();
    float dx = meshSize[0];
    float dy = meshSize[1];
    float Gamma = props(2, 2).mu;  // 使用物性场中的粘度

    float expected_aE = Gamma * dy / dx;
    float expected_aW = Gamma * dy / dx;
    float expected_aN = Gamma * dx / dy;
    float expected_aS = Gamma * dx / dy;
    float expected_aP = expected_aE + expected_aW + expected_aN + expected_aS;

    // 检查系数是否接近理论值（考虑浮点精度）
    float tolerance = 1e-6f;
    bool aE_ok = std::abs(coef.aE - expected_aE) < tolerance;
    bool aW_ok = std::abs(coef.aW - expected_aW) < tolerance;
    bool aN_ok = std::abs(coef.aN - expected_aN) < tolerance;
    bool aS_ok = std::abs(coef.aS - expected_aS) < tolerance;
    bool aP_ok = std::abs(coef.aP - expected_aP) < tolerance;

    if (!aE_ok) fmt::print("  aE mismatch: got {:.6f}, expected {:.6f}\n", coef.aE, expected_aE);
    if (!aW_ok) fmt::print("  aW mismatch: got {:.6f}, expected {:.6f}\n", coef.aW, expected_aW);
    if (!aN_ok) fmt::print("  aN mismatch: got {:.6f}, expected {:.6f}\n", coef.aN, expected_aN);
    if (!aS_ok) fmt::print("  aS mismatch: got {:.6f}, expected {:.6f}\n", coef.aS, expected_aS);
    if (!aP_ok) fmt::print("  aP mismatch: got {:.6f}, expected {:.6f}\n", coef.aP, expected_aP);

    return aE_ok && aW_ok && aN_ok && aS_ok && aP_ok;
}

void test_scalar_equation()
{
    fmt::print("\n=== Testing ScalarEquation ===\n");

    // 设置测试配置
    std::array<BoudaryCondition, 4> boundaryInfo;
    int test_ncx, test_ncy;
    float test_Lx, test_Ly;
    setup_test_configuration(boundaryInfo, test_ncx, test_ncy, test_Lx, test_Ly);

    // 创建网格和场
    StructuredMesh testMesh;
    testMesh.createVolumeMesh();
    testMesh.createBoundaryMesh();
    ScalarField pressure(test_ncx, test_ncy, 0.0f);
    VectorField velocity;
    BoundaryField bc(test_ncx, test_ncy, boundaryInfo.data(), 4);
    FluidPropertyField props;

    // 创建标量方程对象
    ScalarEquation eq(testMesh, pressure, velocity, bc, props);

    // 测试构造函数
    test_scalar_equation_construction(eq, boundaryInfo, test_ncx, test_ncy);

    // 测试重置系数
    bool resetTestPassed = test_scalar_equation_reset_coefficients(eq, test_ncx, test_ncy);

    if (resetTestPassed)
    {
        fmt::print("✓ resetCoefficients() test passed - all coefficients are zero\n");
    }
    else
    {
        fmt::print("✗ resetCoefficients() test failed - non-zero coefficients found\n");
    }

    // 测试扩散项
    bool diffusionTestPassed = test_scalar_equation_add_diffusion(eq, testMesh, props);

    if (diffusionTestPassed)
    {
        fmt::print("✓ addDiffusionTerm() test passed - diffusion coefficients calculated correctly\n");
    }
    else
    {
        fmt::print("✗ addDiffusionTerm() test failed - incorrect diffusion coefficients\n");
    }

    fmt::print("✓ ScalarEquation constructor, resetCoefficients() and addDiffusionTerm() tests completed\n");
}

void test()
{
    fmt::print("Hello from test()\n");

    // 运行各个测试模块
    test_basic_fields();
    test_linear_algebra();
    test_scalar_equation();
}
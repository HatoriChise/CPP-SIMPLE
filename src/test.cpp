// src/test.cpp

#include "test.hpp"

void test()
{
    fmt::print("Hello from test()\n");

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




    // 线性代数示例
    Eigen::MatrixXd A = Eigen::MatrixXd::Random(3, 3);
    Eigen::VectorXd b(3);
    b << 1, 2, 3;
    Eigen::VectorXd x = A.colPivHouseholderQr().solve(b);
    fmt::print("A = \n{}\n", A);
    fmt::print("Solution x = [{:.4f}, {:.4f}, {:.4f}]\n", x(0), x(1), x(2));

    double pi = boost::math::constants::pi<double>();
    fmt::print("pi ≈ {:.10f}\n", pi);
}
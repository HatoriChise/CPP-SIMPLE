// src/test.cpp

#include "test.hpp"

void test()
{
    // 使用 fmt 替代 cout
    fmt::print("🌟 Hello from vcpkg-managed dependencies!\n");

    // 使用 Eigen
    Eigen::MatrixXd A = Eigen::MatrixXd::Random(3, 3);
    Eigen::VectorXd b(3);
    b << 1, 2, 3;
    Eigen::VectorXd x = A.colPivHouseholderQr().solve(b);

    fmt::print("A = \n{}\n", A);
    fmt::print("Solution x = [{:.4f}, {:.4f}, {:.4f}]\n", x(0), x(1), x(2));

    // 使用 Boost
    double pi = boost::math::constants::pi<double>();
    fmt::print("π ≈ {:.10f}\n", pi);
}
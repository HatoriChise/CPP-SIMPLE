#pragma once

#include "scalarField.hpp"

#include <string>

struct Velocity
{
    float u;
    float v;
};

// 速度场，由两个标量场组成，不处理边界
class VectorField
{
public:
    VectorField(int ncx, int ncy, const std::string &name = "velocity");

    // 只读访问，(i, j) -> (u_[j][i], v_[j][i])
    Velocity operator()(int i, int j) const;

    ScalarField &u() { return u_; }
    ScalarField &v() { return v_; }
    const ScalarField &u() const { return u_; }
    const ScalarField &v() const { return v_; }
    const std::string &name() const { return name_; }

private:
    std::string name_{};
    ScalarField u_;
    ScalarField v_;
};

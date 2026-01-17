#include "vectorField.hpp"

VectorField::VectorField()
    : u_(::ncx, ::ncy, initalVelocity[0]),
      v_(::ncx, ::ncy, initalVelocity[1])
{
    // 构造时使用 config 中的网格尺寸与初始速度值
}

Velocity VectorField::operator()(int i, int j) const
{
    // (i, j) -> data_[j][i]，与 ScalarField 的映射保持一致
    return Velocity{u_(i, j), v_(i, j)};
}

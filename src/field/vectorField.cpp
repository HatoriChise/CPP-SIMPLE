#include "vectorField.hpp"

VectorField::VectorField(int ncx, int ncy, const std::string &name)
    : name_(name), u_(ncx, ncy, name + "_u"), v_(ncx, ncy, name + "_v")
{
}

Velocity VectorField::operator()(int i, int j) const
{
    // (i, j) -> data_[j][i]，与 ScalarField 的映射保持一致
    return Velocity{u_(i, j), v_(i, j)};
}

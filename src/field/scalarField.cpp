#include "scalarField.hpp"

#include <algorithm>

ScalarField::ScalarField(int ncx, int ncy, const std::string &name)
    : ncx_(ncx), ncy_(ncy), name_(name), data_(boost::extents[ncy][ncx])
{
    assert(ncx_ > 0 && ncy_ > 0);
}

float &ScalarField::operator()(int i, int j)
{
    // (i, j) -> data_[j][i]，j 为行(y)，i 为列(x)，保证内层循环 i 时内存连续
    assert(i >= 0 && i < ncx_);
    assert(j >= 0 && j < ncy_);
    return data_[j][i];
}

const float &ScalarField::operator()(int i, int j) const
{
    assert(i >= 0 && i < ncx_);
    assert(j >= 0 && j < ncy_);
    return data_[j][i];
}

void ScalarField::fill(float value)
{
    std::fill(data_.data(), data_.data() + data_.num_elements(), value);
}

float *ScalarField::data_ptr()
{
    return data_.data();
}

const float *ScalarField::data_ptr() const
{
    return data_.data();
}

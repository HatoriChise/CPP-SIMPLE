#include "scalarField.hpp"

#include <algorithm>

ScalarField::ScalarField() : ncx_(::ncx), ncy_(::ncy), data_(boost::extents[ncy_][ncx_])
{
    assert(ncx_ > 0 && ncy_ > 0);
    fill(initalTemperature); // 默认填充初始温度
}

ScalarField::ScalarField(int ncx, int ncy, float value) : ncx_(ncx), ncy_(ncy), data_(boost::extents[ncy_][ncx_])
{
    assert(ncx_ > 0 && ncy_ > 0);
    fill(value);
}

float &ScalarField::operator()(int i, int j)
{
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

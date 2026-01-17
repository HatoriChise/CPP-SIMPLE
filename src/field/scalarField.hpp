#pragma once

#include <boost/multi_array.hpp>
#include "../config/simulationConfig.hpp"

#include <cassert>
#include <string>

// 标量场，仅负责数据存储，不处理边界
class ScalarField
{
public:
    // 直接使用 config 的网格尺寸与初始温度构造
    ScalarField();

    // 用于速度场等初始值构造
    ScalarField(int ncx, int ncy, float init_value);

    // 读写访问，(i, j) -> data_[j][i]，确保 x 方向连续访问
    float &operator()(int i, int j);
    // 只读访问，(i, j) -> data_[j][i]
    const float &operator()(int i, int j) const;

    void fill(float value);

    float *data_ptr();
    const float *data_ptr() const;

    int ncx() const
    {
        return ncx_;
    }
    int ncy() const
    {
        return ncy_;
    }

private:
    int ncx_{};
    int ncy_{};
    boost::multi_array<float, 2> data_;
};

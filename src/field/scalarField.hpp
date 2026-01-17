#pragma once

#include <boost/multi_array.hpp>

#include <cassert>
#include <string>

// 标量场，仅负责数据存储，不处理边界
class ScalarField
{
public:
    ScalarField(int ncx, int ncy, const std::string &name = "");

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
    const std::string &name() const
    {
        return name_;
    }

private:
    int ncx_{};
    int ncy_{};
    std::string name_{};
    boost::multi_array<float, 2> data_;
};

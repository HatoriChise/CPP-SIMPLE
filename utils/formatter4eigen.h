// utils/formatter4eigen.h
#pragma once
#include <Eigen/Dense>
#include <fmt/ostream.h>

// Enable fmt formatting for all Eigen matrices/vectors via their operator<<
template<typename Scalar, int Rows, int Cols, int MaxRows, int MaxCols>
struct fmt::formatter<Eigen::Matrix<Scalar, Rows, Cols, MaxRows, MaxCols>>
    : fmt::ostream_formatter
{
    /* override */
};

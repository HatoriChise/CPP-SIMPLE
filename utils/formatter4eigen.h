// utils/formatter4eigen.h
#pragma once
#include <Eigen/Dense>
#include <fmt/ostream.h>
#include "../src/field/vectorField.hpp"

// Enable fmt formatting for all Eigen matrices/vectors via their operator<<
template<typename Scalar, int Rows, int Cols, int MaxRows, int MaxCols>
struct fmt::formatter<Eigen::Matrix<Scalar, Rows, Cols, MaxRows, MaxCols>>
    : fmt::ostream_formatter
{
    /* override */
};

template <>
struct fmt::formatter<Velocity> {
    // 解析格式说明符（如 {:.2f}），此处支持基础格式转发
    char presentation = 'f';
    int precision = 4;

    constexpr auto parse(format_parse_context& ctx) -> decltype(ctx.begin()) {
        auto it = ctx.begin(), end = ctx.end();
        if (it != end && (*it == 'f' || *it == 'e' || *it == 'g')) {
            presentation = *it++;
        }
        return it;
    }

    template <typename FormatContext>
    auto format(const Velocity& v, FormatContext& ctx) const -> decltype(ctx.out()) {

        return fmt::format_to(ctx.out(), "[u: {:.4f}, v: {:.4f}]", v.u, v.v);
    }
};
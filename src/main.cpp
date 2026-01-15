// src/main.cpp
#include <fmt/format.h>
#include <chrono>

#include "src/test.hpp"

int main()
{
    auto t0 = std::chrono::high_resolution_clock::now();

    test();

    // end timing
    auto t1 = std::chrono::high_resolution_clock::now();
    using ms = std::chrono::duration<double, std::milli>;
    fmt::print("Elapsed time: {:.2f} ms\n", std::chrono::duration_cast<ms>(t1 - t0).count());
    return 0;
}
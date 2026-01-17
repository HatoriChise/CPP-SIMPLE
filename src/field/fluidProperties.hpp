#pragma once

#include "scalarField.hpp"
#include "../config/simulationConfig.hpp"

#include <cassert>
#include <string>

// 单点物性集合，便于一次性返回
struct FluidValues {
    float rho;
    float mu;
    float k;
    float cp;
    float nu;    // 运动粘度 = mu / rho
    float alpha; // 热扩散率 = k / (rho * cp)
};

// 物性场：每个 cell 存储物性，可扩展到变物性（温度/压力相关）
// 不做边界处理，外部根据需要写入/更新
class FluidPropertyField {
public:
    FluidPropertyField(int ncx, int ncy, const std::string& base_name = "fluid");

    // 一次性填充全域常量物性（仍作为场存储，便于未来变物性）
    void fill(float rho, float mu, float k, float cp);
    // 直接使用 config 中的物性常量进行填充（mu 由 rho*nu 推得）
    void fill_from_config();

    // 访问单元物性，(i, j) -> data[j][i] 映射与 ScalarField 一致
    FluidValues operator()(int i, int j) const;

    // 分量访问，方便外部根据温度场等自行更新
    ScalarField& rho() { return rho_; }
    ScalarField& mu() { return mu_; }
    ScalarField& k() { return k_; }
    ScalarField& cp() { return cp_; }
    const ScalarField& rho() const { return rho_; }
    const ScalarField& mu() const { return mu_; }
    const ScalarField& k() const { return k_; }
    const ScalarField& cp() const { return cp_; }

    int ncx() const { return ncx_; }
    int ncy() const { return ncy_; }
    const std::string& name() const { return name_; }

private:
    int ncx_{};
    int ncy_{};
    std::string name_{};

    ScalarField rho_; // 密度场
    ScalarField mu_;  // 动力粘度场
    ScalarField k_;   // 导热率场
    ScalarField cp_;  // 定压比热场
};

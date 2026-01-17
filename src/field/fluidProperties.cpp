// src/field/fluidProperties.cpp

#include "fluidProperties.hpp"
#include <iostream>

FluidPropertyField::FluidPropertyField(int ncx, int ncy, const std::string& base_name)
    : ncx_(ncx),
      ncy_(ncy),
      name_(base_name),
      rho_(ncx, ncy, base_name + "_rho"),
      mu_(ncx, ncy, base_name + "_mu"),
      k_(ncx, ncy, base_name + "_k"),
      cp_(ncx, ncy, base_name + "_cp")
{
    // 初始化时，默认设为 0 或某个物理上有意义的极小值，防止除零错误
    fill(1.0f, 1e-5f, 0.026f, 1005.0f); 
}

void FluidPropertyField::fill(float rho_val, float mu_val, float k_val, float cp_val)
{
    rho_.fill(rho_val);
    mu_.fill(mu_val);
    k_.fill(k_val);
    cp_.fill(cp_val);
}

void FluidPropertyField::fill_from_config()
{
    // 从 simulationConfig.hpp 获取常量
    float mu_val = nu * density; 
    
    // 假设 config 中也有导热率 k 和比热 cp，如果没有，可暂设默认值
    // 这里使用你 config 中定义的常量
    fill(density, mu_val, k_fluid, cp_fluid);
    
    std::cout << "FluidPropertyField initialized from config: "
              << "rho=" << density << ", mu=" << mu_val << std::endl;
}

FluidValues FluidPropertyField::operator()(int i, int j) const
{
    // 这里 operator() 内部会调用 ScalarField 的 operator()，映射为 [j][i]
    float current_rho = rho_(i, j);
    float current_mu  = mu_(i, j);
    float current_k   = k_(i, j);
    float current_cp  = cp_(i, j);

    FluidValues values;
    values.rho   = current_rho;
    values.mu    = current_mu;
    values.k     = current_k;
    values.cp    = current_cp;
    
    // 计算导出量
    values.nu    = current_mu / current_rho;
    values.alpha = current_k / (current_rho * current_cp);

    return values;
}
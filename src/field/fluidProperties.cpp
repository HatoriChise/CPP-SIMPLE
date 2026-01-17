// src/field/fluidProperties.cpp

#include "fluidProperties.hpp"
#include <iostream>

FluidPropertyField::FluidPropertyField()
    : ncx_(::ncx),
      ncy_(::ncy),
      rho_(ncx_, ncy_, density),
      mu_(ncx_, ncy_, mu_fluid),
      k_(ncx_, ncy_, k_fluid),
      cp_(ncx_, ncy_, cp_fluid)
{
    // 初始化时，默认设为 0 或某个物理上有意义的极小值，防止除零错误
}

void FluidPropertyField::fill(float rho_val, float mu_val, float k_val, float cp_val)
{
    rho_.fill(rho_val);
    mu_.fill(mu_val);
    k_.fill(k_val);
    cp_.fill(cp_val);
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
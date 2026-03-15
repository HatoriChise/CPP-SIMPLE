// src/euation/scalarEquation.cpp

#include "src/equation/scalarEquation.hpp"
#include "src/field/massFluxField.hpp"

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/SparseCholesky>

#include <fmt/core.h>
#include <boost/math/constants/constants.hpp>

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <vector>

// Eigen 格式化输出
#include "utils/formatter4eigen.h"

// ============================================================================
// 构造与析构
// ============================================================================

ScalarEquation::ScalarEquation(StructuredMesh &mesh, ScalarField &scalarField,
                               VectorField &vectorField, BoundaryField &boundaryField,
                               FluidPropertyField &fluidPropertyField, int direction)
    : mesh_(mesh), scalarField_(scalarField), vectorField_(vectorField),
      boundaryField_(boundaryField), fluidPropertyField_(fluidPropertyField), direction_(direction)
{
    // 初始化系数矩阵
    coefMatrix_.resize(boost::extents[ncy][ncx]);
    resetCoefficients();
}

ScalarEquation::~ScalarEquation()
{
    // 无需特殊清理
}

// ============================================================================
// 系数矩阵操作
// ============================================================================

void ScalarEquation::applyBoundaryCondition()
{
    // direction_ = 0 是 u动量方程，direction_ = 1 是 v动量方程，-1 是标量方程
    // coefficients : aP * phi_P = aE * phi_E + aW * phi_W + aN * phi_N + aS * phi_S + bsrc
    // aE aW aN aS

    // get mesh size
    auto meshSize = mesh_.getMeshSize();
    float dx = meshSize[0];
    float dy = meshSize[1];

    // 遍历所有边界单元，根据边界条件调整系数
    // west i = 0
    for(int j = 0; j < ncy; ++j)
    {
        const auto &bc = boundaryField_.west(j);
        float &aW = coefMatrix_[j][0].aW;
        float &aP = coefMatrix_[j][0].aP;
        float &bsrc = coefMatrix_[j][0].bsrc;

        if(bc.velocityType == INLET)
        {
            // 入口：速度已知，使用 Dirichlet 条件
            float u_inlet = bc.VelocityValue[0];
            float v_inlet = bc.VelocityValue[1];
            // 低阶精度插值
            if(direction_ == 0)
            {
                aP += aW;
                bsrc += 2.0f * aW * u_inlet;
                aW = 0.0f;
            }
            if (direction_ == 1)
            {
                aP += aW;
                bsrc += 2.0f * aW * v_inlet;
                aW = 0.0f;
            }
        }
        else if(bc.velocityType == OUTLET)
        {
            // 出口：速度未知，使用 Neumann 条件（零梯度）
            aP -= aW;
            aW = 0.0f;
        }
        else if(bc.velocityType == WALL)
        {
            // 壁面：速度为零，使用 Dirichlet 条件
            aP += aW;
            aW = 0.0f;
        }
    }

    // East i = ncx - 1
    for(int j = 0; j < ncy; ++j)
    {
        const auto &bc = boundaryField_.east(j);
        float &aE = coefMatrix_[j][ncx - 1].aE;
        float &aP = coefMatrix_[j][ncx - 1].aP;
        float &bsrc = coefMatrix_[j][ncx - 1].bsrc;

        if(bc.velocityType == INLET)
        {
            // 入口：速度已知，使用 Dirichlet 条件
            float u_inlet = bc.VelocityValue[0];
            float v_inlet = bc.VelocityValue[1];
            if(direction_ == 0)
            {
                aP += aE;
                bsrc += 2.0f * aE * u_inlet;
                aE = 0.0f;
            }
            if (direction_ == 1)
            {
                aP += aE;
                bsrc += 2.0f * aE * v_inlet;
                aE = 0.0f;
            }
        }
        else if(bc.velocityType == OUTLET)
        {
            // 出口：速度未知，使用 Neumann 条件（零梯度）
            aP -= aE;
            aE = 0.0f;
        }
        else if(bc.velocityType == WALL)
        {
            // 壁面：速度为零，使用 Dirichlet 条件
            aP += aE;
            aE = 0.0f;
        }
    }

    // South j = 0
    for(int i = 0; i < ncx; ++i)
    {
        const auto &bc = boundaryField_.south(i);
        float &aS = coefMatrix_[0][i].aS;
        float &aP = coefMatrix_[0][i].aP;
        float &bsrc = coefMatrix_[0][i].bsrc;

        if(bc.velocityType == INLET)
        {
            // 入口：速度已知，使用 Dirichlet 条件
            float u_inlet = bc.VelocityValue[0];
            float v_inlet = bc.VelocityValue[1];
            if(direction_ == 0)
            {
                aP += aS;
                bsrc += 2.0f * aS * u_inlet;
                aS = 0.0f;
            }
            if (direction_ == 1)
            {
                aP += aS;
                bsrc += 2.0f * aS * v_inlet;
                aS = 0.0f;
            }
        }
        else if(bc.velocityType == OUTLET)
        {
            // 出口：速度未知，使用 Neumann 条件（零梯度）
            aP -= aS;
            aS = 0.0f;
        }
        else if(bc.velocityType == WALL)
        {
            // 壁面：速度为零，使用 Dirichlet 条件
            aP += aS;
            aS = 0.0f;
        }
    }

    // North j = ncy - 1
    for(int i = 0; i < ncx; ++i)
    {
        const auto &bc = boundaryField_.north(i);
        float &aN = coefMatrix_[ncy - 1][i].aN;
        float &aP = coefMatrix_[ncy - 1][i].aP;
        float &bsrc = coefMatrix_[ncy - 1][i].bsrc;

        if(bc.velocityType == INLET)
        {
            // 入口：速度已知，使用 Dirichlet 条件
            float u_inlet = bc.VelocityValue[0];
            float v_inlet = bc.VelocityValue[1];
            if(direction_ == 0)
            {
                aP += aN;
                bsrc += 2.0f * aN * u_inlet;
                aN = 0.0f;
            }
            if (direction_ == 1)
            {
                aP += aN;
                bsrc += 2.0f * aN * v_inlet;
                aN = 0.0f;
            }
        }
        else if(bc.velocityType == OUTLET)
        {
            // 出口：速度未知，使用 Neumann 条件（零梯度）
            aP -= aN;
            aN = 0.0f;
        }
        else if(bc.velocityType == WALL)
        {
            // 壁面：速度为零，使用 Dirichlet 条件
            aP += aN;
            aN = 0.0f;
        }
    }
}

void ScalarEquation::resetCoefficients()
{
    // 遍历所有单元
    for(int j = 0; j < ncy; ++j)
    {
        for(int i = 0; i < ncx; ++i)
        {
            coefMatrix_[j][i].aE = 0.0f;
            coefMatrix_[j][i].aW = 0.0f;
            coefMatrix_[j][i].aN = 0.0f;
            coefMatrix_[j][i].aS = 0.0f;
            coefMatrix_[j][i].aP = 0.0f;
            coefMatrix_[j][i].bsrc = 0.0f;
        }
    }
}

// ============================================================================
// 界面通量计算
// ============================================================================

// 计算界面扩散系数（调和平均）
float ScalarEquation::computeFaceDiffusionCoefficient(float mu_owner, float mu_neighbor,
                                                      float distance, float area) const
{
    // 调和平均：mu_face = 2 * mu_P * mu_N / (mu_P + mu_N)
    // 如果 mu_P + mu_N 接近 0，避免除零
    float mu_sum = mu_owner + mu_neighbor;
    if(mu_sum < 1e-10f)
    {
        return 0.0f;
    }

    float mu_face = 2.0f * mu_owner * mu_neighbor / mu_sum;
    return mu_face * area / distance;
}

void ScalarEquation::addDiffusionTerm()
{
    // 获取单元中心坐标和网格尺寸
    const auto &xc = mesh_.getCellCentersX();
    const auto &yc = mesh_.getCellCentersY();
    auto meshSize = mesh_.getMeshSize();
    float dx = meshSize[0]; // x 方向单元尺寸（用于北/南面面积）
    float dy = meshSize[1]; // y 方向单元尺寸（用于东/西面面积）

    // 遍历所有单元（包括边界单元）
    for(int j = 0; j < ncy; ++j)
    {
        for(int i = 0; i < ncx; ++i)
        {
            float mu_P = fluidPropertyField_(i, j).mu;

            // 计算相邻单元中心距离
            float dx_PE = (i + 1 < ncx) ? (xc[i + 1] - xc[i]) : dx; // P 到 E 的距离
            float dx_WP = (i - 1 >= 0) ? (xc[i] - xc[i - 1]) : dx;  // W 到 P 的距离
            float dy_PN = (j + 1 < ncy) ? (yc[j + 1] - yc[j]) : dy; // P 到 N 的距离
            float dy_SP = (j - 1 >= 0) ? (yc[j] - yc[j - 1]) : dy;  // S 到 P 的距离

            // 东界面 (e)
            float aE = 0.0f;
            if(i + 1 < ncx)
            {
                // 内部界面：使用调和平均
                float mu_E = fluidPropertyField_(i + 1, j).mu;
                aE = computeFaceDiffusionCoefficient(mu_P, mu_E, dx_PE, dy);
            }
            else
            {
                // 东边界：虚拟单元 mu_E = mu_P
                aE = computeFaceDiffusionCoefficient(mu_P, mu_P, dx_PE, dy);
            }

            // 西界面 (w)
            float aW = 0.0f;
            if(i - 1 >= 0)
            {
                // 内部界面：使用调和平均
                float mu_W = fluidPropertyField_(i - 1, j).mu;
                aW = computeFaceDiffusionCoefficient(mu_P, mu_W, dx_WP, dy);
            }
            else
            {
                // 西边界：虚拟单元 mu_W = mu_P
                aW = computeFaceDiffusionCoefficient(mu_P, mu_P, dx_WP, dy);
            }

            // 北界面 (n)
            float aN = 0.0f;
            if(j + 1 < ncy)
            {
                // 内部界面：使用调和平均
                float mu_N = fluidPropertyField_(i, j + 1).mu;
                aN = computeFaceDiffusionCoefficient(mu_P, mu_N, dy_PN, dx);
            }
            else
            {
                // 北边界：虚拟单元 mu_N = mu_P
                aN = computeFaceDiffusionCoefficient(mu_P, mu_P, dy_PN, dx);
            }

            // 南界面 (s)
            float aS = 0.0f;
            if(j - 1 >= 0)
            {
                // 内部界面：使用调和平均
                float mu_S = fluidPropertyField_(i, j - 1).mu;
                aS = computeFaceDiffusionCoefficient(mu_P, mu_S, dy_SP, dx);
            }
            else
            {
                // 南边界：虚拟单元 mu_S = mu_P
                aS = computeFaceDiffusionCoefficient(mu_P, mu_P, dy_SP, dx);
            }

            // 累加到系数矩阵
            coefMatrix_[j][i].aE += aE;
            coefMatrix_[j][i].aW += aW;
            coefMatrix_[j][i].aN += aN;
            coefMatrix_[j][i].aS += aS;
            coefMatrix_[j][i].aP += (aE + aW + aN + aS);
        }
    }
}


void ScalarEquation::addConvectionTerm(const MassFluxField &massFlux)
{
    // 遍历所有单元（包括边界）
    for(int j = 0; j < ncy; ++j)
    {
        for(int i = 0; i < ncx; ++i)
        {
            // 直接由存储场读取已知通量：东面和北面直接读取
            float F_e = massFlux(i, j).mE;
            float F_n = massFlux(i, j).mN;
            float F_w = massFlux(i, j).mW;
            float F_s = massFlux(i, j).mS;

            // 迎风格式计算对流系数
            coefMatrix_[j][i].aE += std::max(-F_e, 0.0f);
            coefMatrix_[j][i].aW += std::max(F_w, 0.0f);
            coefMatrix_[j][i].aN += std::max(-F_n, 0.0f);
            coefMatrix_[j][i].aS += std::max(F_s, 0.0f);

            // 中心系数贡献（连续性）
            coefMatrix_[j][i].aP += std::max(-F_e, 0.0f) + std::max(F_w, 0.0f) +
                                    std::max(-F_n, 0.0f) + std::max(F_s, 0.0f);
        }
    }
}

void ScalarEquation::addSourceTerm()
{
    // TODO: 实现源项
}

void ScalarEquation::addPressureGradient(const ScalarField &pressure)
{
    
    if (direction_ < 0)
    {
        fmt::print("Warning: addPressureGradient called for scalar equation. No action taken.\n");
        return;
    }

    // 获取网格尺寸
    auto meshSize = mesh_.getMeshSize();
    float dx = meshSize[0];
    float dy = meshSize[1];

    float area_x = dy; // x 方向界面面积
    float area_y = dx; // y 方向界面面积
    float volyme = dx * dy; // 单元体积

    const auto& bc = boundaryField_;

    for (int j = 0; j < ncy; ++j)
    {
        for (int i = 0; i < ncx; ++i)
        {
            // x direction pressure gradient
            if (direction_ == 0)
            {
                float pl, pr; // pl = pressure at west face, pr = pressure at east face

                const auto& bc_west = bc.west(j);
                const auto& bc_east = bc.east(j);

                // West face
                if (i == 0) 
                {
                    // 西边界：根据边界条件计算 pl
                    if (bc_west.velocityType == INLET || bc_west.velocityType == WALL)
                    {
                        pl = pressure(i, j); // 虚拟单元压力等于当前单元压力
                    }
                    else if (bc_west.velocityType == OUTLET)
                    {
                        pl = bc_west.pressureValue; // 使用指定压力值
                    }
                    else
                    {
                        pl = pressure(i, j); // 默认处理
                    }
                }
                else
                {
                    pl = pressure(i - 1, j); // 内部单元压力
                }

                if (i = ncx - 1)
                {
                    // 东边界：根据边界条件计算 pr
                    if (bc_east.velocityType == INLET || bc_east.velocityType == WALL)
                    {
                        // 壁面/入口：零梯度 ∂p/∂x = 0 → p[ncx] = p[ncx-1]
                        pr = pressure(i, j); // 虚拟单元压力等于当前单元压力
                    }
                    else if (bc_east.velocityType == OUTLET)
                    {
                        pr = bc_east.pressureValue; // 使用指定压力值
                    }
                    else
                    {
                        pr = pressure(i, j); // 默认处理
                    }
                }
                else
                {
                    pr = pressure(i + 1, j); // 内部单元压力
                }

                // 添加压力梯度源项：S = - (pr - pl) / (2 * dx) * volume
                coefMatrix_[j][i].bsrc += 0.5f * (pr - pl) * area_x;
            }

            if (direction_ == 1)
            {
                float pl, pr; // pl = p[j-1] at south face, pr = p[j+1] at north face

                const auto& bc_south = bc.south(i);
                const auto& bc_north = bc.north(i);

                // South face
                if (j == 0) 
                {
                    // 南边界：根据边界条件计算 pb
                    if (bc_south.velocityType == INLET || bc_south.velocityType == WALL)
                    {
                        pl = pressure(i, j); // 虚拟单元压力等于当前单元压力
                    }
                    else if (bc_south.velocityType == OUTLET)
                    {
                        pl = bc_south.pressureValue; // 使用指定压力值
                    }
                    else
                    {
                        pl = pressure(i, j); // 默认处理
                    }
                }
                else
                {
                    pl = pressure(i, j - 1); // 内部单元压力
                }

                if (j == ncy - 1)
                {
                    // 北边界：根据边界条件计算 pt
                    if (bc_north.velocityType == INLET || bc_north.velocityType == WALL)
                    {
                        pr = pressure(i, j); // 虚拟单元压力等于当前单元压力
                    }
                    else if (bc_north.velocityType == OUTLET)
                    {
                        pr = bc_north.pressureValue; // 使用指定压力值
                    }
                    else
                    {
                        pr = pressure(i, j); // 默认处理
                    }
                }
                else
                {
                    pr = pressure(i, j + 1); // 内部单元压力
                }

                // 添加压力梯度源项：S = - (pr - pl) / (2 * dy) * volume
                coefMatrix_[j][i].bsrc += 0.5f * (pr - pl) * area_y;
            }
        }
    }
}

void ScalarEquation::setRelaxation(float relaxationFactor)
{
    // 亚松驰：a_P = a_P / relaxationFactor
    // 源项修正：b = b + (1-relaxationFactor)*a_P*phi_P / relaxationFactor
    for(int j = 0; j < ncy; ++j)
    {
        for(int i = 0; i < ncx; ++i)
        {
            float aP_old = coefMatrix_[j][i].aP;
            coefMatrix_[j][i].aP /= relaxationFactor;
            // 源项修正（如果需要考虑亚松驰）
            // coefMatrix_[j][i].bsrc += (1.0f - relaxationFactor) * aP_old * scalarField_(i, j) /
            // relaxationFactor;
        }
    }
}

// ============================================================================
// 工具函数
// ============================================================================

void ScalarEquation::saveCoefficientsToFile(const std::string &filename) const
{
    std::ofstream file(filename);
    if(!file.is_open())
    {
        std::cerr << "Failed to open file: " << filename << std::endl;
        return;
    }

    // CSV 头
    file << "i,j,aE,aW,aN,aS,aP,bsrc\n";

    // 写入系数
    for(int j = 0; j < ncy; ++j)
    {
        for(int i = 0; i < ncx; ++i)
        {
            const auto &coef = coefMatrix_[j][i];
            file << i << "," << j << "," << coef.aE << "," << coef.aW << "," << coef.aN << ","
                 << coef.aS << "," << coef.aP << "," << coef.bsrc << "\n";
        }
    }

    file.close();
    fmt::print("Coefficients saved to {}\n", filename);
}

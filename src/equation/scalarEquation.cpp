// src/euation/scalarEquation.cpp

#include "src/equation/scalarEquation.hpp"

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

bool ScalarEquation::isBoundaryFace(int i, int j, Face face) const
{
    switch(face)
    {
    case Face::East: // 东面 (e)
        return (i + 1 >= ncx);
    case Face::West: // 西面 (w)
        return (i - 1 < 0);
    case Face::North: // 北面 (n)
        return (j + 1 >= ncy);
    case Face::South: // 南面 (s)
        return (j - 1 < 0);
    default:
        throw std::invalid_argument("Invalid face index");
    }
}

float ScalarEquation::computeFaceMassFlux(int i, int j, Face face,
                                          const ScalarField *pressure) const
{
    auto meshSize = mesh_.getMeshSize();
    float dx = meshSize[0]; // x 方向单元尺寸（用于北/南面面积）
    float dy = meshSize[1]; // y 方向单元尺寸（用于东/西面面积）
    float volume = dx * dy; // 单元体积
    float rho = fluidPropertyField_(i, j).rho;


    float linearVelocity = 0.0f;
    float area = (face == Face::East || face == Face::West) ? dy : dx;
    switch(face)
    {
    case Face::East: // 东面 (e)
        linearVelocity = 0.5f * (vectorField_.u()(i, j) + vectorField_.u()(i + 1, j));
        break;
    case Face::West: // 西面 (w)
        linearVelocity = 0.5f * (vectorField_.u()(i - 1, j) + vectorField_.u()(i, j));
        break;
    case Face::North: // 北面 (n)
        linearVelocity = 0.5f * (vectorField_.v()(i, j) + vectorField_.v()(i, j + 1));
        break;
    case Face::South: // 南面 (s)
        linearVelocity = 0.5f * (vectorField_.v()(i, j - 1) + vectorField_.v()(i, j));
        break;
    default:
        throw std::invalid_argument("Invalid face index");
    }

    if(pressure == nullptr)
    {
        // 线性插值
        return rho * linearVelocity * area;
    }

    // RC interpolation
    // 计算压力梯度修正项
    float aP = coefMatrix_[j][i].aP;

    float aP_nabor = 0.0f;
    switch(face)
    {
    case Face::East: // 东面 (e)
        aP_nabor = coefMatrix_[j][i + 1].aP;
        break;
    case Face::West: // 西面 (w)
        aP_nabor = coefMatrix_[j][i - 1].aP;
        break;
    case Face::North: // 北面 (n)
        aP_nabor = coefMatrix_[j + 1][i].aP;
        break;
    case Face::South: // 南面 (s)
        aP_nabor = coefMatrix_[j - 1][i].aP;
        break;
    }

    // 防止 aP 过小导致数值不稳定，回退到线性插值通量
    if(std::fabs(aP) < 1e-10f || std::fabs(aP_nabor) < 1e-10f)
    {
        return rho * linearVelocity * area;
    }

    // 计算系数d = V / aP
    float d_P = volume / aP;
    float d_nabor = volume / aP_nabor;
    float d_face = 0.5f * (d_P + d_nabor);

    // RC插值公式
    /**
     * $$u_e=\bar{u}_e+\frac{1}{2}d_P\left( \frac{P_E-P_W}{2\Delta x} \right) +\frac{1}{2}d_E\left(
     * \frac{P_{EE}-P_P}{2\Delta x} \right) -d_e\left( \frac{P_E - P_P}{\Delta x} \right) $$
     */
    // 计算(P_E - P_W) / (2 * dx)
    float gradP_P = 0.0f;
    switch(face)
    {
    case Face::East: // 东面 (e)
    case Face::West: // 西面 (w)
        if(i > 0 && i < ncx - 1)
        {
            gradP_P = (pressure->operator()(i + 1, j) - pressure->operator()(i - 1, j)) / (2.0f * dx);
        }
        else if(i == 0)
        {
            gradP_P = (pressure->operator()(i + 1, j) - pressure->operator()(i, j)) / dx;
        }
        else
        {
            gradP_P = (pressure->operator()(i, j) - pressure->operator()(i - 1, j)) / dx;
        }
        break;
    case Face::North: // 北面 (n)
    case Face::South: // 南面 (s)
        if(j > 0 && j < ncy - 1)
        {
            gradP_P = (pressure->operator()(i, j + 1) - pressure->operator()(i, j - 1)) / (2.0f * dy);
        }
        else if(j == 0)
        {
            gradP_P = (pressure->operator()(i, j + 1) - pressure->operator()(i, j)) / dy;
        }
        else
        {
            gradP_P = (pressure->operator()(i, j) - pressure->operator()(i, j - 1)) / dy;
        }
        break;
    }
    
    // 计算(P_EE - P_P) / (2 * dx)
    float gradP_nabor = 0.0f;
    switch(face)
    {
    case Face::East: // 邻居为 E
        if(i < ncx - 2)
        {
            // 在 E 邻居点做中心差分：(P_{i+2} - P_i) / (2*dx)
            gradP_nabor = (pressure->operator()(i + 2, j) - pressure->operator()(i, j)) / (2.0f * dx);
        }
        else
        {
            // E 邻居落在右边界附近，回退到单边差分
            gradP_nabor = (pressure->operator()(i + 1, j) - pressure->operator()(i, j)) / dx;
        }
        break;
    case Face::West: // 邻居为 W
        if(i > 1)
        {
            // 在 W 邻居点做中心差分：(P_i - P_{i-2}) / (2*dx)
            gradP_nabor = (pressure->operator()(i, j) - pressure->operator()(i - 2, j)) / (2.0f * dx);
        }
        else
        {
            // W 邻居落在左边界附近，回退到单边差分
            gradP_nabor = (pressure->operator()(i, j) - pressure->operator()(i - 1, j)) / dx;
        }
        break;
    case Face::North: // 邻居为 N
        if(j < ncy - 2)
        {
            // 在 N 邻居点做中心差分：(P_{j+2} - P_j) / (2*dy)
            gradP_nabor = (pressure->operator()(i, j + 2) - pressure->operator()(i, j)) / (2.0f * dy);
        }
        else
        {
            // N 邻居落在上边界附近，回退到单边差分
            gradP_nabor = (pressure->operator()(i, j + 1) - pressure->operator()(i, j)) / dy;
        }
        break;
    case Face::South: // 邻居为 S
        if(j > 1)
        {
            // 在 S 邻居点做中心差分：(P_j - P_{j-2}) / (2*dy)
            gradP_nabor = (pressure->operator()(i, j) - pressure->operator()(i, j - 2)) / (2.0f * dy);
        }
        else
        {
            // S 邻居落在下边界附近，回退到单边差分
            gradP_nabor = (pressure->operator()(i, j) - pressure->operator()(i, j - 1)) / dy;
        }
        break;
    }
    
    // 计算(P_E - P_P) / dx
    float gradP_face = 0.0f;
    switch(face)
    {
    case Face::East:
        if(i + 1 < ncx)
        {
            gradP_face = (pressure->operator()(i + 1, j) - pressure->operator()(i, j)) / dx;
        }
        else if(i - 1 >= 0)
        {
            gradP_face = (pressure->operator()(i, j) - pressure->operator()(i - 1, j)) / dx;
        }
        else
        {
            gradP_face = 0.0f;
        }
        break;
    case Face::West:
        if(i - 1 >= 0)
        {
            gradP_face = (pressure->operator()(i, j) - pressure->operator()(i - 1, j)) / dx;
        }
        else if(i + 1 < ncx)
        {
            gradP_face = (pressure->operator()(i + 1, j) - pressure->operator()(i, j)) / dx;
        }
        else
        {
            gradP_face = 0.0f;
        }
        break;
    case Face::North:
        if(j + 1 < ncy)
        {
            gradP_face = (pressure->operator()(i, j + 1) - pressure->operator()(i, j)) / dy;
        }
        else if(j - 1 >= 0)
        {
            gradP_face = (pressure->operator()(i, j) - pressure->operator()(i, j - 1)) / dy;
        }
        else
        {
            gradP_face = 0.0f;
        }
        break;
    case Face::South:
        if(j - 1 >= 0)
        {
            gradP_face = (pressure->operator()(i, j) - pressure->operator()(i, j - 1)) / dy;
        }
        else if(j + 1 < ncy)
        {
            gradP_face = (pressure->operator()(i, j + 1) - pressure->operator()(i, j)) / dy;
        }
        else
        {
            gradP_face = 0.0f;
        }
        break;
    default:
        throw std::invalid_argument("Invalid face index");
    }

    float correctedVelocity =
        linearVelocity + 0.5f * d_P * gradP_P + 0.5f * d_nabor * gradP_nabor - d_face * gradP_face;

    return rho * correctedVelocity * area;
}

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

void ScalarEquation::addConvectionTerm(const ScalarField *pressure)
{
    // 遍历所有单元（包括边界）
    // computeFaceMassFlux 内部会处理边界情况（返回 0 通量）
    for(int j = 0; j < ncy; ++j)
    {
        for(int i = 0; i < ncx; ++i)
        {
            // 计算四个界面的质量通量
            float F_e = computeFaceMassFlux(i, j, Face::East, pressure);  // 东
            float F_w = computeFaceMassFlux(i, j, Face::West, pressure);  // 西
            float F_n = computeFaceMassFlux(i, j, Face::North, pressure); // 北
            float F_s = computeFaceMassFlux(i, j, Face::South, pressure); // 南

            // 迎风格式计算对流系数
            coefMatrix_[j][i].aE += std::max(-F_e, 0.0f);
            coefMatrix_[j][i].aW += std::max(F_w, 0.0f);
            coefMatrix_[j][i].aN += std::max(-F_n, 0.0f);
            coefMatrix_[j][i].aS += std::max(F_s, 0.0f);

            // 中心系数贡献（连续性）
            coefMatrix_[j][i].aP += (F_e - F_w + F_n - F_s);
        }
    }
}

void ScalarEquation::addSourceTerm()
{
    // TODO: 实现源项
}

void ScalarEquation::addPressureGradient(const ScalarField &pressure)
{
    // 获取网格尺寸
    auto meshSize = mesh_.getMeshSize();
    float dx = meshSize[0];
    float dy = meshSize[1];

    // 遍历内部单元
    for(int j = 1; j < ncy - 1; ++j)
    {
        for(int i = 1; i < ncx - 1; ++i)
        {
            if(direction_ == 0)
            { // u 动量方程：-dp/dx * V
                float p_e = 0.5f * (pressure(i, j) + pressure(i + 1, j));
                float p_w = 0.5f * (pressure(i - 1, j) + pressure(i, j));
                float dpdx = (p_e - p_w) / dx;
                float volume = dx * dy;
                coefMatrix_[j][i].bsrc -= dpdx * volume;
            }
            else if(direction_ == 1)
            { // v 动量方程：-dp/dy * V
                float p_n = 0.5f * (pressure(i, j) + pressure(i, j + 1));
                float p_s = 0.5f * (pressure(i, j - 1) + pressure(i, j));
                float dpdy = (p_n - p_s) / dy;
                float volume = dx * dy;
                coefMatrix_[j][i].bsrc -= dpdy * volume;
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

void ScalarEquation::applyBoundaries()
{
    // TODO: 实现边界条件应用
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

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


void ScalarEquation::addConvectionTerm(const MassFluxField& massFlux)
{
    // 遍历所有单元（包括边界）
    for(int j = 0; j < ncy; ++j)
    {
        for(int i = 0; i < ncx; ++i)
        {
            // 直接由存储场读取已知通量：东面和北面直接读取
            float F_e = massFlux(i, j).mE;
            float F_n = massFlux(i, j).mN;
            
            // 西面和南面优先从外部专门的通量库读取
            float F_w = massFlux(i, j).mW;
            float F_s = massFlux(i, j).mS;

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

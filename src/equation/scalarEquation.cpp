// src/euation/scalarEquation.cpp

#include "src/equation/scalarEquation.hpp"

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/SparseCholesky>

#include <boost/math/constants/constants.hpp>
#include <fmt/core.h>

#include <algorithm>
#include <fstream>
#include <iostream>
#include <limits>
#include <vector>

// Eigen 格式化输出
#include "utils/formatter4eigen.h"

// 网格尺寸（从配置中获取）
extern int ncx, ncy;
extern float Lx, Ly;

// ============================================================================
// 构造与析构
// ============================================================================

ScalarEquation::ScalarEquation(StructuredMesh &mesh, ScalarField &scalarField,
                               VectorField &vectorField, BoundaryField &boundaryField,
                               FluidPropertyField &fluidPropertyField, int direction)
    : mesh_(mesh), scalarField_(scalarField), vectorField_(vectorField),
      boundaryField_(boundaryField), fluidPropertyField_(fluidPropertyField),
      direction_(direction)
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
    for (int j = 0; j < ncy; ++j)
    {
        for (int i = 0; i < ncx; ++i)
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

float ScalarEquation::computeFaceMassFlux(int i, int j, int face, const ScalarField* pressure) const
{
    auto meshSize = mesh_.getMeshSize();
    float dx = meshSize[0];
    float dy = meshSize[1];
    float rho = fluidPropertyField_(i, j).rho;

    float u_face = 0.0f, v_face = 0.0f, area = 0.0f;
    float dp = 0.0f;  // 压力差（用于 Rhie-Chow 修正）

    switch (face)
    {
    case 0:  // 东侧 (i+1/2, j)
        u_face = 0.5f * (vectorField_.u()(i, j) + vectorField_.u()(i + 1, j));
        v_face = 0.5f * (vectorField_.v()(i, j) + vectorField_.v()(i + 1, j));
        area = dy;
        if (pressure) dp = (*pressure)(i, j) - (*pressure)(i + 1, j);  // p_P - p_E
        break;
    case 1:  // 西侧 (i-1/2, j)
        u_face = 0.5f * (vectorField_.u()(i - 1, j) + vectorField_.u()(i, j));
        v_face = 0.5f * (vectorField_.v()(i - 1, j) + vectorField_.v()(i, j));
        area = dy;
        if (pressure) dp = (*pressure)(i - 1, j) - (*pressure)(i, j);  // p_W - p_P
        break;
    case 2:  // 北侧 (i, j+1/2)
        u_face = 0.5f * (vectorField_.u()(i, j) + vectorField_.u()(i, j + 1));
        v_face = 0.5f * (vectorField_.v()(i, j) + vectorField_.v()(i, j + 1));
        area = dx;
        if (pressure) dp = (*pressure)(i, j) - (*pressure)(i, j + 1);  // p_P - p_N
        break;
    case 3:  // 南侧 (i, j-1/2)
        u_face = 0.5f * (vectorField_.u()(i, j - 1) + vectorField_.u()(i, j));
        v_face = 0.5f * (vectorField_.v()(i, j - 1) + vectorField_.v()(i, j));
        area = dx;
        if (pressure) dp = (*pressure)(i, j - 1) - (*pressure)(i, j);  // p_S - p_P
        break;
    }

    // Rhie-Chow 修正（如果提供了压力场且 aP>0）
    if (pressure != nullptr && coefMatrix_[j][i].aP > 0.0f)
    {
        float d = area / coefMatrix_[j][i].aP;  // 动量系数倒数
        // 根据界面方向修正对应速度分量
        if (face == 0 || face == 1)  // 东/西界面，修正 u
        {
            u_face += d * dp;
        }
        else  // 北/南界面，修正 v
        {
            v_face += d * dp;
        }
    }

    // 返回质量通量
    if (face == 0 || face == 1)
        return rho * u_face * area;
    else
        return rho * v_face * area;
}

// 计算界面扩散系数（调和平均）
float ScalarEquation::computeFaceDiffusionCoefficient(float mu_owner, float mu_neighbor, 
                                                       float distance, float area) const
{
    // 调和平均：mu_face = 2 * mu_P * mu_N / (mu_P + mu_N)
    // 如果 mu_P + mu_N 接近 0，避免除零
    float mu_sum = mu_owner + mu_neighbor;
    if (mu_sum < 1e-10f)
    {
        return 0.0f;
    }
    
    float mu_face = 2.0f * mu_owner * mu_neighbor / mu_sum;
    return mu_face * area / distance;
}

void ScalarEquation::addDiffusionTerm()
{
    // 获取单元中心坐标和网格尺寸
    const auto& xc = mesh_.getCellCentersX();
    const auto& yc = mesh_.getCellCentersY();
    auto meshSize = mesh_.getMeshSize();
    float dx = meshSize[0];  // x 方向单元尺寸（用于北/南面面积）
    float dy = meshSize[1];  // y 方向单元尺寸（用于东/西面面积）

    // 遍历所有单元（包括边界单元）
    for (int j = 0; j < ncy; ++j)
    {
        for (int i = 0; i < ncx; ++i)
        {
            float mu_P = fluidPropertyField_(i, j).mu;

            // 计算相邻单元中心距离
            float dx_PE = (i + 1 < ncx) ? (xc[i + 1] - xc[i]) : dx;  // P 到 E 的距离
            float dx_WP = (i - 1 >= 0) ? (xc[i] - xc[i - 1]) : dx;   // W 到 P 的距离
            float dy_PN = (j + 1 < ncy) ? (yc[j + 1] - yc[j]) : dy;  // P 到 N 的距离
            float dy_SP = (j - 1 >= 0) ? (yc[j] - yc[j - 1]) : dy;   // S 到 P 的距离

            // 东界面 (e)
            float aE = 0.0f;
            if (i + 1 < ncx)
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
            if (i - 1 >= 0)
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
            if (j + 1 < ncy)
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
            if (j - 1 >= 0)
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

void ScalarEquation::addConvectionTerm(const ScalarField* pressure)
{
    // 遍历内部单元
    for (int j = 1; j < ncy - 1; ++j)
    {
        for (int i = 1; i < ncx - 1; ++i)
        {
            // 计算界面质量通量（传递 pressure 以启用 Rhie-Chow 插值）
            float F_e = computeFaceMassFlux(i, j, 0, pressure);  // 东 East
            float F_w = computeFaceMassFlux(i, j, 1, pressure);  // 西 West
            float F_n = computeFaceMassFlux(i, j, 2, pressure);  // 北 North
            float F_s = computeFaceMassFlux(i, j, 3, pressure);  // 南 South

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

void ScalarEquation::addPressureGradient(const ScalarField& pressure)
{
    // 获取网格尺寸
    auto meshSize = mesh_.getMeshSize();
    float dx = meshSize[0];
    float dy = meshSize[1];

    // 遍历内部单元
    for (int j = 1; j < ncy - 1; ++j)
    {
        for (int i = 1; i < ncx - 1; ++i)
        {
            if (direction_ == 0)
            {  // u 动量方程：-dp/dx * V
                float p_e = 0.5f * (pressure(i, j) + pressure(i + 1, j));
                float p_w = 0.5f * (pressure(i - 1, j) + pressure(i, j));
                float dpdx = (p_e - p_w) / dx;
                float volume = dx * dy;
                coefMatrix_[j][i].bsrc -= dpdx * volume;
            }
            else if (direction_ == 1)
            {  // v 动量方程：-dp/dy * V
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
    for (int j = 0; j < ncy; ++j)
    {
        for (int i = 0; i < ncx; ++i)
        {
            float aP_old = coefMatrix_[j][i].aP;
            coefMatrix_[j][i].aP /= relaxationFactor;
            // 源项修正（如果需要考虑亚松驰）
            // coefMatrix_[j][i].bsrc += (1.0f - relaxationFactor) * aP_old * scalarField_(i, j) / relaxationFactor;
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

void ScalarEquation::saveCoefficientsToFile(const std::string& filename) const
{
    std::ofstream file(filename);
    if (!file.is_open())
    {
        std::cerr << "Failed to open file: " << filename << std::endl;
        return;
    }

    // CSV 头
    file << "i,j,aE,aW,aN,aS,aP,bsrc\n";

    // 写入系数
    for (int j = 0; j < ncy; ++j)
    {
        for (int i = 0; i < ncx; ++i)
        {
            const auto& coef = coefMatrix_[j][i];
            file << i << "," << j << "," << coef.aE << "," << coef.aW << ","
                 << coef.aN << "," << coef.aS << "," << coef.aP << "," << coef.bsrc << "\n";
        }
    }

    file.close();
    fmt::print("Coefficients saved to {}\n", filename);
}

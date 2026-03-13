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

    float u_face = 0.0f, v_face = 0.0f;
    float area = 0.0f;

    // 线性插值速度
    float u_linear = 0.0f, v_linear = 0.0f;

    switch (face)
    {
    case 0:  // 东侧 (i+1/2, j)
    {
        u_linear = 0.5f * (vectorField_.u()(i, j) + vectorField_.u()(i + 1, j));
        v_linear = 0.5f * (vectorField_.v()(i, j) + vectorField_.v()(i + 1, j));
        area = dy;

        // RC 修正条件：有压力场、当前单元 aP 有效
        bool apply_rc = (pressure != nullptr && coefMatrix_[j][i].aP > 0.0f);

        if (apply_rc)
        {
            // 获取 d 系数：d = area / aP
            float d_P = area / coefMatrix_[j][i].aP;
            
            // 检查东侧邻居是否存在
            bool has_E = (i + 1 < ncx);
            float d_E = has_E ? (area / coefMatrix_[j][i + 1].aP) : d_P;

            // 计算压力梯度
            // gradP_P: P 单元中心压力梯度
            float gradP_P = 0.0f;
            bool has_W = (i - 1 >= 0);
            if (has_E && has_W)
            {
                // 中心差分：(P_E - P_W) / (2*dx)
                gradP_P = ((*pressure)(i + 1, j) - (*pressure)(i - 1, j)) / (2.0f * dx);
            }
            else if (has_E && !has_W)
            {
                // 西边界：向前差分 (P_E - P_P) / dx
                gradP_P = ((*pressure)(i + 1, j) - (*pressure)(i, j)) / dx;
            }
            else if (!has_E && has_W)
            {
                // 东边界：向后差分 (P_P - P_W) / dx
                gradP_P = ((*pressure)(i, j) - (*pressure)(i - 1, j)) / dx;
            }
            else
            {
                // 孤立单元（不应发生）
                gradP_P = 0.0f;
            }

            // gradP_E: E 单元中心压力梯度
            float gradP_E = 0.0f;
            bool has_EE = has_E && (i + 2 < ncx);
            if (has_E)
            {
                if (has_EE)
                {
                    // 中心差分：(P_EE - P_P) / (2*dx)
                    gradP_E = ((*pressure)(i + 2, j) - (*pressure)(i, j)) / (2.0f * dx);
                }
                else
                {
                    // 东边界附近：向前差分 (P_E - P_P) / dx
                    gradP_E = ((*pressure)(i + 1, j) - (*pressure)(i, j)) / dx;
                }
            }

            // gradP_e: 界面压力梯度
            float gradP_e = 0.0f;
            if (has_E)
            {
                gradP_e = ((*pressure)(i, j) - (*pressure)(i + 1, j)) / dx;
            }
            else
            {
                // 东边界：使用向后差分
                gradP_e = ((*pressure)(i, j) - (*pressure)(i - 1, j)) / dx;
            }

            // Rhie-Chow 修正
            // u_e = ū_e + 0.5*d_P*gradP_P + 0.5*d_E*gradP_E - d_e*gradP_e
            // 其中 d_e = 0.5*(d_P + d_E)
            float d_e = 0.5f * (d_P + d_E);
            u_face = u_linear + 0.5f * d_P * gradP_P + 0.5f * d_E * gradP_E - d_e * gradP_e;
            v_face = v_linear;
        }
        else
        {
            u_face = u_linear;
            v_face = v_linear;
        }
        break;
    }

    case 1:  // 西侧 (i-1/2, j)
    {
        u_linear = 0.5f * (vectorField_.u()(i - 1, j) + vectorField_.u()(i, j));
        v_linear = 0.5f * (vectorField_.v()(i - 1, j) + vectorField_.v()(i, j));
        area = dy;

        // RC 修正条件：有压力场、当前单元和西侧单元 aP 都有效
        bool apply_rc = (pressure != nullptr &&
                         coefMatrix_[j][i].aP > 0.0f &&
                         (i - 1 >= 0) && coefMatrix_[j][i - 1].aP > 0.0f);

        if (apply_rc)
        {
            // 获取 W 和 P 单元的 aP 系数
            float aP_W = coefMatrix_[j][i - 1].aP;
            float aP_P = coefMatrix_[j][i].aP;

            // 计算 d 系数：d = area / aP
            float d_W = area / aP_W;
            float d_P = area / aP_P;
            float d_w = 0.5f * (d_W + d_P);

            // 检查西侧是否有邻居 (W-1)
            bool has_WW = (i - 2 >= 0);

            // 计算单元中心压力梯度 (中心差分)
            float dpdx_W = 0.0f;
            if (has_WW)
            {
                // (∂p/∂x)_W = (p_P - p_{W-1}) / (2*dx)
                dpdx_W = (pressure->operator()(i, j) - pressure->operator()(i - 2, j)) / (2.0f * dx);
            }
            else
            {
                // 西边界：单侧差分 (p_P - p_W) / dx
                dpdx_W = (pressure->operator()(i, j) - pressure->operator()(i - 1, j)) / dx;
            }

            // (∂p/∂x)_P = (p_E - p_W) / (2*dx)
            bool has_E = (i + 1 < ncx);
            float dpdx_P = 0.0f;
            if (has_E)
            {
                dpdx_P = (pressure->operator()(i + 1, j) - pressure->operator()(i - 1, j)) / (2.0f * dx);
            }
            else
            {
                // 东边界：单侧差分 (p_P - p_W) / dx
                dpdx_P = (pressure->operator()(i, j) - pressure->operator()(i - 1, j)) / dx;
            }

            // 界面压力梯度 (∂p/∂x)_w = (p_W - p_P) / dx
            float dpdx_w = (pressure->operator()(i - 1, j) - pressure->operator()(i, j)) / dx;

            // Rhie-Chow 修正
            u_face = u_linear + 0.5f * d_W * dpdx_W + 0.5f * d_P * dpdx_P - d_w * dpdx_w;
            v_face = v_linear;
        }
        else
        {
            u_face = u_linear;
            v_face = v_linear;
        }
        break;
    }

    case 2:  // 北侧 (i, j+1/2)
    {
        u_linear = 0.5f * (vectorField_.u()(i, j) + vectorField_.u()(i, j + 1));
        v_linear = 0.5f * (vectorField_.v()(i, j) + vectorField_.v()(i, j + 1));
        area = dx;

        // RC 修正条件：有压力场、当前单元 aP 有效
        bool apply_rc = (pressure != nullptr && coefMatrix_[j][i].aP > 0.0f);

        if (apply_rc)
        {
            float d_P = area / coefMatrix_[j][i].aP;
            
            bool has_N = (j + 1 < ncy);
            float d_N = has_N ? (area / coefMatrix_[j + 1][i].aP) : d_P;

            // 计算 y 方向压力梯度
            float gradP_P = 0.0f;
            bool has_S = (j - 1 >= 0);
            if (has_N && has_S)
            {
                gradP_P = ((*pressure)(i, j + 1) - (*pressure)(i, j - 1)) / (2.0f * dy);
            }
            else if (has_N && !has_S)
            {
                gradP_P = ((*pressure)(i, j + 1) - (*pressure)(i, j)) / dy;
            }
            else if (!has_N && has_S)
            {
                gradP_P = ((*pressure)(i, j) - (*pressure)(i, j - 1)) / dy;
            }

            float gradP_N = 0.0f;
            bool has_NN = has_N && (j + 2 < ncy);
            if (has_N)
            {
                if (has_NN)
                {
                    gradP_N = ((*pressure)(i, j + 2) - (*pressure)(i, j)) / (2.0f * dy);
                }
                else
                {
                    gradP_N = ((*pressure)(i, j + 1) - (*pressure)(i, j)) / dy;
                }
            }

            float gradP_n = 0.0f;
            if (has_N)
            {
                gradP_n = ((*pressure)(i, j) - (*pressure)(i, j + 1)) / dy;
            }
            else
            {
                gradP_n = ((*pressure)(i, j) - (*pressure)(i, j - 1)) / dy;
            }

            float d_n = 0.5f * (d_P + d_N);
            v_face = v_linear + 0.5f * d_P * gradP_P + 0.5f * d_N * gradP_N - d_n * gradP_n;
            u_face = u_linear;
        }
        else
        {
            u_face = u_linear;
            v_face = v_linear;
        }
        break;
    }

    case 3:  // 南侧 (i, j-1/2)
    {
        u_linear = 0.5f * (vectorField_.u()(i, j - 1) + vectorField_.u()(i, j));
        v_linear = 0.5f * (vectorField_.v()(i, j - 1) + vectorField_.v()(i, j));
        area = dx;

        bool apply_rc = (pressure != nullptr &&
                         coefMatrix_[j][i].aP > 0.0f &&
                         (j - 1 >= 0) && coefMatrix_[j][j - 1].aP > 0.0f);

        if (apply_rc)
        {
            float aP_S = coefMatrix_[j][j - 1].aP;
            float aP_P = coefMatrix_[j][i].aP;

            float d_S = area / aP_S;
            float d_P = area / aP_P;
            float d_s = 0.5f * (d_S + d_P);

            // 检查南侧是否有邻居 (S-1)
            bool has_SS = (j - 2 >= 0);
            float dpdy_S = 0.0f;
            if (has_SS)
            {
                dpdy_S = (pressure->operator()(i, j) - pressure->operator()(i, j - 2)) / (2.0f * dy);
            }
            else
            {
                // 南边界：单侧差分
                dpdy_S = (pressure->operator()(i, j) - pressure->operator()(i, j - 1)) / dy;
            }

            bool has_N = (j + 1 < ncy);
            float dpdy_P = 0.0f;
            if (has_N)
            {
                dpdy_P = (pressure->operator()(i, j + 1) - pressure->operator()(i, j - 1)) / (2.0f * dy);
            }
            else
            {
                // 北边界：单侧差分
                dpdy_P = (pressure->operator()(i, j) - pressure->operator()(i, j - 1)) / dy;
            }

            float dpdy_s = (pressure->operator()(i, j - 1) - pressure->operator()(i, j)) / dy;

            v_face = v_linear + 0.5f * d_S * dpdy_S + 0.5f * d_P * dpdy_P - d_s * dpdy_s;
            u_face = u_linear;
        }
        else
        {
            u_face = u_linear;
            v_face = v_linear;
        }
        break;
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
    // 遍历所有单元（包括边界）
    for (int j = 0; j < ncy; ++j)
    {
        for (int i = 0; i < ncx; ++i)
        {
            // 计算界面质量通量（传递 pressure 以启用 Rhie-Chow 插值）
            // 边界处需要检查邻居是否存在
            float F_e = 0.0f, F_w = 0.0f, F_n = 0.0f, F_s = 0.0f;

            // 东界面：需要东侧有邻居
            if (i + 1 < ncx)
            {
                F_e = computeFaceMassFlux(i, j, 0, pressure);
            }

            // 西界面：需要西侧有邻居
            if (i - 1 >= 0)
            {
                F_w = computeFaceMassFlux(i, j, 1, pressure);
            }

            // 北界面：需要北侧有邻居
            if (j + 1 < ncy)
            {
                F_n = computeFaceMassFlux(i, j, 2, pressure);
            }

            // 南界面：需要南侧有邻居
            if (j - 1 >= 0)
            {
                F_s = computeFaceMassFlux(i, j, 3, pressure);
            }

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

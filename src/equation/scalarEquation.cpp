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

    // 检查界面是否为域边界（边界处无通量）
    bool is_domain_boundary = false;

    switch (face)
    {
    case 0:  // 东侧 (i+1/2, j)
        // 如果东侧是域边界，返回 0 通量
        if (i + 1 >= ncx - 1)
        {
            is_domain_boundary = true;
        }
        else
        {
            u_linear = 0.5f * (vectorField_.u()(i, j) + vectorField_.u()(i + 1, j));
            v_linear = 0.5f * (vectorField_.v()(i, j) + vectorField_.v()(i + 1, j));
            area = dy;
        }
        break;

    case 1:  // 西侧 (i-1/2, j)
        // 如果西侧是域边界，返回 0 通量
        if (i - 1 <= 0)
        {
            is_domain_boundary = true;
        }
        else
        {
            u_linear = 0.5f * (vectorField_.u()(i - 1, j) + vectorField_.u()(i, j));
            v_linear = 0.5f * (vectorField_.v()(i - 1, j) + vectorField_.v()(i, j));
            area = dy;
        }
        break;

    case 2:  // 北侧 (i, j+1/2)
        // 如果北侧是域边界，返回 0 通量
        if (j + 1 >= ncy - 1)
        {
            is_domain_boundary = true;
        }
        else
        {
            u_linear = 0.5f * (vectorField_.u()(i, j) + vectorField_.u()(i, j + 1));
            v_linear = 0.5f * (vectorField_.v()(i, j) + vectorField_.v()(i, j + 1));
            area = dx;
        }
        break;

    case 3:  // 南侧 (i, j-1/2)
        // 如果南侧是域边界，返回 0 通量
        if (j - 1 <= 0)
        {
            is_domain_boundary = true;
        }
        else
        {
            u_linear = 0.5f * (vectorField_.u()(i, j - 1) + vectorField_.u()(i, j));
            v_linear = 0.5f * (vectorField_.v()(i, j - 1) + vectorField_.v()(i, j));
            area = dx;
        }
        break;
    }

    // 域边界处无通量
    if (is_domain_boundary)
    {
        return 0.0f;
    }

    // RC 修正（如果提供了压力场）
    // Rhie-Chow 插值公式：u_face = ū_face + 0.5*d_P*(∇p)_P + 0.5*d_N*(∇p)_N - d_face*(∇p)_face
    // 其中：
    //   d_P = area / aP_P, d_N = area / aP_N, d_face = 0.5*(d_P + d_N)
    //   (∇p)_P: P 单元中心压力梯度
    //   (∇p)_N: 邻居单元中心压力梯度（东界面为 E，西界面为 W，北界面为 N，南界面为 S）
    //   (∇p)_face: 界面压力梯度
    if (pressure != nullptr && coefMatrix_[j][i].aP > 0.0f)
    {
        // 根据界面方向获取相邻单元的 aP
        float aP_P = coefMatrix_[j][i].aP;
        float aP_N = 0.0f;  // 邻居单元的 aP

        switch (face)
        {
        case 0:  // 东侧：邻居为 E (i+1, j)
            aP_N = coefMatrix_[j][i + 1].aP;
            break;
        case 1:  // 西侧：邻居为 W (i-1, j)
            aP_N = coefMatrix_[j][i - 1].aP;
            break;
        case 2:  // 北侧：邻居为 N (i, j+1)
            aP_N = coefMatrix_[j + 1][i].aP;
            break;
        case 3:  // 南侧：邻居为 S (i, j-1)
            aP_N = coefMatrix_[j - 1][i].aP;
            break;
        }

        // 如果邻居 aP 无效，不使用 RC
        if (aP_N <= 0.0f)
        {
            if (face == 0 || face == 1)
                return rho * u_linear * area;
            else
                return rho * v_linear * area;
        }

        // 计算 d 系数：d = area / aP
        float d_P = area / aP_P;
        float d_N = area / aP_N;  // 东侧为 d_E，西侧为 d_W，北侧为 d_N，南侧为 d_S
        float d_face = 0.5f * (d_P + d_N);

        // 计算压力梯度
        float gradP_P = 0.0f;   // P 单元中心压力梯度
        float gradP_N = 0.0f;   // 邻居单元中心压力梯度
        float gradP_face = 0.0f;  // 界面压力梯度

        switch (face)
        {
        case 0:  // 东侧 (i+1/2, j) - x 方向
        {
            // === Rhie-Chow 插值公式（东界面）===
            // u_e = ū_e + 0.5*d_P*(∇p)_P + 0.5*d_E*(∇p)_E - d_e*(∇p)_e
            // 其中：
            //   (∇p)_P ≈ (P_E - P_W) / (2Δx)   [P 单元中心压力梯度，中心差分]
            //   (∇p)_E ≈ (P_EE - P_P) / (2Δx)  [E 单元中心压力梯度，中心差分]
            //   (∇p)_e ≈ (P_E - P_P) / Δx      [界面 e 压力梯度]
            
            // 检查邻居是否存在
            bool has_W = (i - 1 >= 0);
            bool has_E = (i + 1 < ncx - 1);
            bool has_EE = (i + 2 < ncx - 1);

            // 计算 (∇p)_P：P 单元中心压力梯度
            if (has_W && has_E)
            {
                // 中心差分：(P_E - P_W) / (2Δx)
                gradP_P = ((*pressure)(i + 1, j) - (*pressure)(i - 1, j)) / (2.0f * dx);
            }
            else if (has_E && !has_W)
            {
                // 西边界附近：向前差分 (P_E - P_P) / Δx
                gradP_P = ((*pressure)(i + 1, j) - (*pressure)(i, j)) / dx;
            }
            else if (has_W && !has_E)
            {
                // 东边界附近：向后差分 (P_P - P_W) / Δx
                gradP_P = ((*pressure)(i, j) - (*pressure)(i - 1, j)) / dx;
            }

            // 计算 (∇p)_E：E 单元中心压力梯度（注意：这是 gradP_E，不是 gradP_N）
            if (has_EE)
            {
                // 中心差分：(P_EE - P_P) / (2Δx)
                gradP_N = ((*pressure)(i + 2, j) - (*pressure)(i, j)) / (2.0f * dx);
            }
            else if (has_E)
            {
                // 东边界附近：向前差分 (P_E - P_P) / Δx
                gradP_N = ((*pressure)(i + 1, j) - (*pressure)(i, j)) / dx;
            }

            // 计算 (∇p)_e：界面 e 压力梯度
            // 公式：(∇p)_e ≈ (P_E - P_P) / Δx
            gradP_face = ((*pressure)(i + 1, j) - (*pressure)(i, j)) / dx;

            // RC 修正
            u_linear += 0.5f * d_P * gradP_P + 0.5f * d_N * gradP_N - d_face * gradP_face;
            break;
        }

        case 1:  // 西侧 (i-1/2, j) - x 方向
        {
            // === Rhie-Chow 插值公式（西界面）===
            // u_w = ū_w + 0.5*d_P*(∇p)_P + 0.5*d_W*(∇p)_W - d_w*(∇p)_w
            // 其中：
            //   (∇p)_P ≈ (P_E - P_W) / (2Δx)   [P 单元中心压力梯度]
            //   (∇p)_W ≈ (P_P - P_WW) / (2Δx)  [W 单元中心压力梯度]
            //   (∇p)_w ≈ (P_W - P_P) / Δx      [界面 w 压力梯度]
            
            bool has_W = (i - 1 >= 0);
            bool has_WW = (i - 2 >= 0);
            bool has_E = (i + 1 < ncx - 1);

            // 计算 (∇p)_P：P 单元中心压力梯度
            if (has_W && has_E)
            {
                gradP_P = ((*pressure)(i + 1, j) - (*pressure)(i - 1, j)) / (2.0f * dx);
            }
            else if (has_W && !has_E)
            {
                gradP_P = ((*pressure)(i, j) - (*pressure)(i - 1, j)) / dx;
            }

            // 计算 (∇p)_W：W 单元中心压力梯度（注意：这是 gradP_W）
            if (has_WW)
            {
                // 中心差分：(P_P - P_WW) / (2Δx)
                gradP_N = ((*pressure)(i, j) - (*pressure)(i - 2, j)) / (2.0f * dx);
            }
            else if (has_W)
            {
                // 西边界附近：向前差分 (P_W - P_P) / Δx（注意符号）
                gradP_N = ((*pressure)(i - 1, j) - (*pressure)(i, j)) / dx;
            }

            // 计算 (∇p)_w：界面 w 压力梯度
            // 公式：(∇p)_w ≈ (P_W - P_P) / Δx
            gradP_face = ((*pressure)(i - 1, j) - (*pressure)(i, j)) / dx;

            // RC 修正
            u_linear += 0.5f * d_P * gradP_P + 0.5f * d_N * gradP_N - d_face * gradP_face;
            break;
        }

        case 2:  // 北侧 (i, j+1/2) - y 方向
        {
            // === Rhie-Chow 插值公式（北界面）===
            // v_n = ū_n + 0.5*d_P*(∇p)_P + 0.5*d_N*(∇p)_N - d_n*(∇p)_n
            // 其中：
            //   (∇p)_P ≈ (P_N - P_S) / (2Δy)   [P 单元中心压力梯度]
            //   (∇p)_N ≈ (P_NN - P_P) / (2Δy)  [N 单元中心压力梯度]
            //   (∇p)_n ≈ (P_N - P_P) / Δy      [界面 n 压力梯度]
            
            bool has_S = (j - 1 >= 0);
            bool has_N = (j + 1 < ncy - 1);
            bool has_NN = (j + 2 < ncy - 1);

            // 计算 (∇p)_P：P 单元中心压力梯度（y 方向）
            if (has_S && has_N)
            {
                gradP_P = ((*pressure)(i, j + 1) - (*pressure)(i, j - 1)) / (2.0f * dy);
            }
            else if (has_N && !has_S)
            {
                gradP_P = ((*pressure)(i, j + 1) - (*pressure)(i, j)) / dy;
            }
            else if (has_S && !has_N)
            {
                gradP_P = ((*pressure)(i, j) - (*pressure)(i, j - 1)) / dy;
            }

            // 计算 (∇p)_N：N 单元中心压力梯度
            if (has_NN)
            {
                gradP_N = ((*pressure)(i, j + 2) - (*pressure)(i, j)) / (2.0f * dy);
            }
            else if (has_N)
            {
                gradP_N = ((*pressure)(i, j + 1) - (*pressure)(i, j)) / dy;
            }

            // 计算 (∇p)_n：界面 n 压力梯度
            // 公式：(∇p)_n ≈ (P_N - P_P) / Δy
            gradP_face = ((*pressure)(i, j + 1) - (*pressure)(i, j)) / dy;

            // RC 修正
            v_linear += 0.5f * d_P * gradP_P + 0.5f * d_N * gradP_N - d_face * gradP_face;
            break;
        }

        case 3:  // 南侧 (i, j-1/2) - y 方向
        {
            // === Rhie-Chow 插值公式（南界面）===
            // v_s = ū_s + 0.5*d_P*(∇p)_P + 0.5*d_S*(∇p)_S - d_s*(∇p)_s
            // 其中：
            //   (∇p)_P ≈ (P_N - P_S) / (2Δy)   [P 单元中心压力梯度]
            //   (∇p)_S ≈ (P_P - P_SS) / (2Δy)  [S 单元中心压力梯度]
            //   (∇p)_s ≈ (P_S - P_P) / Δy      [界面 s 压力梯度]
            
            bool has_S = (j - 1 >= 0);
            bool has_SS = (j - 2 >= 0);
            bool has_N = (j + 1 < ncy - 1);

            // 计算 (∇p)_P：P 单元中心压力梯度（y 方向）
            if (has_S && has_N)
            {
                gradP_P = ((*pressure)(i, j + 1) - (*pressure)(i, j - 1)) / (2.0f * dy);
            }
            else if (has_S && !has_N)
            {
                gradP_P = ((*pressure)(i, j) - (*pressure)(i, j - 1)) / dy;
            }

            // 计算 (∇p)_S：S 单元中心压力梯度
            if (has_SS)
            {
                gradP_N = ((*pressure)(i, j) - (*pressure)(i, j - 2)) / (2.0f * dy);
            }
            else if (has_S)
            {
                gradP_N = ((*pressure)(i, j - 1) - (*pressure)(i, j)) / dy;
            }

            // 计算 (∇p)_s：界面 s 压力梯度
            // 公式：(∇p)_s ≈ (P_S - P_P) / Δy
            gradP_face = ((*pressure)(i, j - 1) - (*pressure)(i, j)) / dy;

            // RC 修正
            v_linear += 0.5f * d_P * gradP_P + 0.5f * d_N * gradP_N - d_face * gradP_face;
            break;
        }
        }
    }

    // 返回质量通量
    if (face == 0 || face == 1)
        return rho * u_linear * area;
    else
        return rho * v_linear * area;
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
    // computeFaceMassFlux 内部会处理边界情况（返回 0 通量）
    for (int j = 0; j < ncy; ++j)
    {
        for (int i = 0; i < ncx; ++i)
        {
            // 计算四个界面的质量通量
            float F_e = computeFaceMassFlux(i, j, 0, pressure);  // 东
            float F_w = computeFaceMassFlux(i, j, 1, pressure);  // 西
            float F_n = computeFaceMassFlux(i, j, 2, pressure);  // 北
            float F_s = computeFaceMassFlux(i, j, 3, pressure);  // 南

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

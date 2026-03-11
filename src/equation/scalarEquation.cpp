#include "scalarEquation.hpp"
#include <fstream>

ScalarEquation::ScalarEquation(StructuredMesh &mesh, ScalarField &scalarField,
                               VectorField &vectorField, BoundaryField &boundaryField,
                               FluidPropertyField &fluidPropertyField, int direction)
    : mesh_(mesh),
      scalarField_(scalarField),
      vectorField_(vectorField),
      boundaryField_(boundaryField),
      fluidPropertyField_(fluidPropertyField),
      direction_(direction)
{
    // 根据全局网格尺寸 (ncx, ncy) 初始化 coefMatrix_ 的大小
    coefMatrix_.resize(boost::extents[ncy][ncx]);

    // 执行一次 resetCoefficients() 确保初始状态为零
    resetCoefficients();
}

ScalarEquation::~ScalarEquation() = default;

void ScalarEquation::resetCoefficients()
{
    // 遍历 coefMatrix_，将所有单元的 aE, aW, aN, aS, aP, bsrc 重置为 0
    // 这是每一轮 SIMPLE 迭代开始前必须执行的操作
    for (int j = 0; j < ncy; ++j)
    {
        for (int i = 0; i < ncx; ++i)
        {
            coefMatrix_[j][i] = COEF{0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f};
        }
    }
}

// === 核心组装方法 ===

float ScalarEquation::computeFaceMassFlux(int i, int j, int face, const ScalarField* pressure) const
{
    auto meshSize = mesh_.getMeshSize();
    float dx = meshSize[0];
    float dy = meshSize[1];
    float volume = dx * dy;  // 单元体积（二维为面积）
    float rho = fluidPropertyField_(i, j).rho;

    float u_face = 0.0f, v_face = 0.0f;
    float area = 0.0f;

    // 线性插值部分的速度
    float u_linear = 0.0f, v_linear = 0.0f;

    // 边界保护标志：在边界附近关闭 RC 修正（相邻单元 aP 可能为零）
    bool near_boundary = (i <= 1 || i >= ncx - 2 || j <= 1 || j >= ncy - 2);

    switch (face)
    {
    case 0:  // 东侧 (i+1/2, j)
    {
        u_linear = 0.5f * (vectorField_.u()(i, j) + vectorField_.u()(i + 1, j));
        v_linear = 0.5f * (vectorField_.v()(i, j) + vectorField_.v()(i + 1, j));
        area = dy;

        // RC 修正条件：有压力场、当前单元和东侧单元 aP 都有效、不在边界附近
        bool apply_rc = (pressure != nullptr && 
                         coefMatrix_[j][i].aP > 0.0f && 
                         coefMatrix_[j][i + 1].aP > 0.0f &&
                         !near_boundary);

        if (apply_rc)
        {
            // 获取 P 和 E 单元的 aP 系数
            float aP_P = coefMatrix_[j][i].aP;
            float aP_E = coefMatrix_[j][i + 1].aP;

            // 计算 d 系数：d = V / aP
            float d_P = volume / aP_P;
            float d_E = volume / aP_E;
            float d_e = 0.5f * (d_P + d_E);

            // 计算单元中心压力梯度 (中心差分)
            // (∂p/∂x)_P = (p_E - p_W) / (2*dx)
            float dpdx_P = (pressure->operator()(i + 1, j) - pressure->operator()(i - 1, j)) / (2.0f * dx);
            // (∂p/∂x)_E = (p_{E+1} - p_P) / (2*dx)
            float dpdx_E = (pressure->operator()(i + 2, j) - pressure->operator()(i, j)) / (2.0f * dx);

            // 界面压力梯度 (∂p/∂x)_e = (p_P - p_E) / dx
            float dpdx_e = (pressure->operator()(i, j) - pressure->operator()(i + 1, j)) / dx;

            // Rhie-Chow 修正后的界面速度
            // u_e = 0.5*(u_P + u_E) + 0.5*d_P*(∂p/∂x)_P + 0.5*d_E*(∂p/∂x)_E - d_e*(∂p/∂x)_e
            u_face = u_linear + 0.5f * d_P * dpdx_P + 0.5f * d_E * dpdx_E - d_e * dpdx_e;
            v_face = v_linear;  // v 分量不需要修正（垂直于界面）
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

        // RC 修正条件：有压力场、当前单元和西侧单元 aP 都有效、不在边界附近
        bool apply_rc = (pressure != nullptr && 
                         coefMatrix_[j][i].aP > 0.0f && 
                         coefMatrix_[j][i - 1].aP > 0.0f &&
                         !near_boundary);

        if (apply_rc)
        {
            // 获取 W 和 P 单元的 aP 系数
            float aP_W = coefMatrix_[j][i - 1].aP;
            float aP_P = coefMatrix_[j][i].aP;

            // 计算 d 系数
            float d_W = volume / aP_W;
            float d_P = volume / aP_P;
            float d_w = 0.5f * (d_W + d_P);

            // 计算单元中心压力梯度 (中心差分)
            // (∂p/∂x)_W = (p_P - p_{W-1}) / (2*dx)
            float dpdx_W = (pressure->operator()(i, j) - pressure->operator()(i - 2, j)) / (2.0f * dx);
            // (∂p/∂x)_P = (p_E - p_W) / (2*dx)
            float dpdx_P = (pressure->operator()(i + 1, j) - pressure->operator()(i - 1, j)) / (2.0f * dx);

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

        // RC 修正条件：有压力场、当前单元和北侧单元 aP 都有效、不在边界附近
        bool apply_rc = (pressure != nullptr && 
                         coefMatrix_[j][i].aP > 0.0f && 
                         coefMatrix_[j + 1][i].aP > 0.0f &&
                         !near_boundary);

        if (apply_rc)
        {
            // 获取 P 和 N 单元的 aP 系数
            float aP_P = coefMatrix_[j][i].aP;
            float aP_N = coefMatrix_[j + 1][i].aP;

            // 计算 d 系数
            float d_P = volume / aP_P;
            float d_N = volume / aP_N;
            float d_n = 0.5f * (d_P + d_N);

            // 计算单元中心压力梯度 (y 方向，中心差分)
            // (∂p/∂y)_P = (p_N - p_S) / (2*dy)
            float dpdy_P = (pressure->operator()(i, j + 1) - pressure->operator()(i, j - 1)) / (2.0f * dy);
            // (∂p/∂y)_N = (p_{N+1} - p_P) / (2*dy)
            float dpdy_N = (pressure->operator()(i, j + 2) - pressure->operator()(i, j)) / (2.0f * dy);

            // 界面压力梯度 (∂p/∂y)_n = (p_P - p_N) / dy
            float dpdy_n = (pressure->operator()(i, j) - pressure->operator()(i, j + 1)) / dy;

            // Rhie-Chow 修正（v 分量）
            v_face = v_linear + 0.5f * d_P * dpdy_P + 0.5f * d_N * dpdy_N - d_n * dpdy_n;
            u_face = u_linear;  // u 分量不需要修正（垂直于界面）
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

        // RC 修正条件：有压力场、当前单元和南侧单元 aP 都有效、不在边界附近
        bool apply_rc = (pressure != nullptr && 
                         coefMatrix_[j][i].aP > 0.0f && 
                         coefMatrix_[j - 1][i].aP > 0.0f &&
                         !near_boundary);

        if (apply_rc)
        {
            // 获取 S 和 P 单元的 aP 系数
            float aP_S = coefMatrix_[j - 1][i].aP;
            float aP_P = coefMatrix_[j][i].aP;

            // 计算 d 系数
            float d_S = volume / aP_S;
            float d_P = volume / aP_P;
            float d_s = 0.5f * (d_S + d_P);

            // 计算单元中心压力梯度 (y 方向，中心差分)
            // (∂p/∂y)_S = (p_P - p_{S-1}) / (2*dy)
            float dpdy_S = (pressure->operator()(i, j) - pressure->operator()(i, j - 2)) / (2.0f * dy);
            // (∂p/∂y)_P = (p_N - p_S) / (2*dy)
            float dpdy_P = (pressure->operator()(i, j + 1) - pressure->operator()(i, j - 1)) / (2.0f * dy);

            // 界面压力梯度 (∂p/∂y)_s = (p_S - p_P) / dy
            float dpdy_s = (pressure->operator()(i, j - 1) - pressure->operator()(i, j)) / dy;

            // Rhie-Chow 修正（v 分量）
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

    default:
        u_face = 0.0f;
        v_face = 0.0f;
        area = 0.0f;
        break;
    }

    // 返回质量通量
    if (face == 0 || face == 1)
        return rho * u_face * area;
    else
        return rho * v_face * area;
}

void ScalarEquation::addDiffusionTerm()
{
    // 获取单元中心坐标和网格尺寸
    const auto& xc = mesh_.getCellCentersX();
    const auto& yc = mesh_.getCellCentersY();
    auto meshSize = mesh_.getMeshSize();
    float dx = meshSize[0];  // x方向单元尺寸（用于北/南面面积）
    float dy = meshSize[1];  // y方向单元尺寸（用于东/西面面积）


    // 遍历内部单元
    for (int j = 1; j < ncy - 1; ++j)
    {
        for (int i = 1; i < ncx - 1; ++i)
        {
            float mu = fluidPropertyField_(i, j).mu;

            // 计算相邻单元中心距离
            float dx_PE = xc[i + 1] - xc[i];  // P到E的距离
            float dx_WP = xc[i] - xc[i - 1];  // W到P的距离
            float dy_PN = yc[j + 1] - yc[j];  // P到N的距离
            float dy_SP = yc[j] - yc[j - 1];  // S到P的距离

            // 扩散系数: a = mu * (面积 / 距离)
            float aE = mu * dy / dx_PE;
            float aW = mu * dy / dx_WP;
            float aN = mu * dx / dy_PN;
            float aS = mu * dx / dy_SP;

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
    // TODO: 将自定义源项场 S 的值映射并累加到 bsrc: bsrc += S(i,j) * volume
    // 针对温度标量方程，S 可以是体积热源项 (W/m^3)，需要乘以单元体积转换为总热量贡献
}

void ScalarEquation::saveCoefficientsToFile(const std::string& filename) const
{
    std::ofstream ofs(filename);
    if (!ofs.is_open())
    {
        fmt::print("ERROR: Cannot open file {} for writing\n", filename);
        return;
    }

    // 写入表头
    ofs << "# i,j,aE,aW,aN,aS,aP,bsrc\n";

    // 遍历所有单元
    for (int j = 0; j < ncy; ++j)
    {
        for (int i = 0; i < ncx; ++i)
        {
            const auto& coef = coefMatrix_[j][i];
            ofs << i << "," << j << ","
                << coef.aE << "," << coef.aW << ","
                << coef.aN << "," << coef.aS << ","
                << coef.aP << "," << coef.bsrc << "\n";
        }
    }

    ofs.close();
    fmt::print("Coefficients saved to {}\n", filename);
}

void ScalarEquation::addPressureGradient(const ScalarField& pressure)
{
    // 仅对动量方程有效
    if (direction_ < 0) return;

    auto meshSize = mesh_.getMeshSize();
    float dx = meshSize[0];
    float dy = meshSize[1];

    // 遍历内部单元
    for (int j = 1; j < ncy - 1; ++j)
    {
        for (int i = 1; i < ncx - 1; ++i)
        {
            float p_P = pressure(i, j);

            if (direction_ == 0)  // u动量方程
            {
                // 界面压力线性插值
                float p_e = 0.5f * (p_P + pressure(i + 1, j));
                float p_w = 0.5f * (pressure(i - 1, j) + p_P);
                // 压力梯度源项
                coefMatrix_[j][i].bsrc -= (p_e - p_w) * dy;
            }
            else if (direction_ == 1)  // v动量方程
            {
                // 界面压力线性插值
                float p_n = 0.5f * (p_P + pressure(i, j + 1));
                float p_s = 0.5f * (pressure(i, j - 1) + p_P);
                // 压力梯度源项
                coefMatrix_[j][i].bsrc -= (p_n - p_s) * dx;
            }
        }
    }
}

void ScalarEquation::applyBoundaries()
{
    // 遍历所有边界单元，应用边界条件
    
    // === 西侧边界 (i=0) ===
    for (int j = 0; j < ncy; ++j)
    {
        const auto& bc = boundaryField_.west(j);
        
        // 处理 u 速度边界条件
        if (direction_ == 0)
        {
            if (bc.velocityType == WALL || bc.velocityType == INLET)
            {
                // Dirichlet 条件：u = u_wall
                float u_bc = bc.VelocityValue[0];
                coefMatrix_[j][0].aP = 1.0f;
                coefMatrix_[j][0].aE = 0.0f;
                coefMatrix_[j][0].aW = 0.0f;
                coefMatrix_[j][0].aN = 0.0f;
                coefMatrix_[j][0].aS = 0.0f;
                coefMatrix_[j][0].bsrc = u_bc;
            }
            else if (bc.velocityType == OUTLET)
            {
                // Neumann 条件：∂u/∂x = 0 → aW = 0, aP = aE
                coefMatrix_[j][0].aW = 0.0f;
                coefMatrix_[j][0].aP = coefMatrix_[j][0].aE;
            }
        }
        
        // 处理 v 速度边界条件
        if (direction_ == 1)
        {
            if (bc.velocityType == WALL || bc.velocityType == INLET)
            {
                // Dirichlet 条件：v = v_wall
                float v_bc = bc.VelocityValue[1];
                coefMatrix_[j][0].aP = 1.0f;
                coefMatrix_[j][0].aE = 0.0f;
                coefMatrix_[j][0].aW = 0.0f;
                coefMatrix_[j][0].aN = 0.0f;
                coefMatrix_[j][0].aS = 0.0f;
                coefMatrix_[j][0].bsrc = v_bc;
            }
        }
    }
    
    // === 东侧边界 (i=ncx-1) ===
    for (int j = 0; j < ncy; ++j)
    {
        const auto& bc = boundaryField_.east(j);
        
        // 处理 u 速度边界条件
        if (direction_ == 0)
        {
            if (bc.velocityType == WALL || bc.velocityType == INLET)
            {
                float u_bc = bc.VelocityValue[0];
                coefMatrix_[j][ncx - 1].aP = 1.0f;
                coefMatrix_[j][ncx - 1].aE = 0.0f;
                coefMatrix_[j][ncx - 1].aW = 0.0f;
                coefMatrix_[j][ncx - 1].aN = 0.0f;
                coefMatrix_[j][ncx - 1].aS = 0.0f;
                coefMatrix_[j][ncx - 1].bsrc = u_bc;
            }
            else if (bc.velocityType == OUTLET)
            {
                // Neumann 条件：∂u/∂x = 0
                coefMatrix_[j][ncx - 1].aE = 0.0f;
                coefMatrix_[j][ncx - 1].aP = coefMatrix_[j][ncx - 1].aW;
            }
        }
        
        // 处理 v 速度边界条件
        if (direction_ == 1)
        {
            if (bc.velocityType == WALL || bc.velocityType == INLET)
            {
                float v_bc = bc.VelocityValue[1];
                coefMatrix_[j][ncx - 1].aP = 1.0f;
                coefMatrix_[j][ncx - 1].aE = 0.0f;
                coefMatrix_[j][ncx - 1].aW = 0.0f;
                coefMatrix_[j][ncx - 1].aN = 0.0f;
                coefMatrix_[j][ncx - 1].aS = 0.0f;
                coefMatrix_[j][ncx - 1].bsrc = v_bc;
            }
        }
    }
    
    // === 南侧边界 (j=0) ===
    for (int i = 0; i < ncx; ++i)
    {
        const auto& bc = boundaryField_.south(i);
        
        // 处理 u 速度边界条件
        if (direction_ == 0)
        {
            if (bc.velocityType == WALL || bc.velocityType == INLET)
            {
                float u_bc = bc.VelocityValue[0];
                coefMatrix_[0][i].aP = 1.0f;
                coefMatrix_[0][i].aE = 0.0f;
                coefMatrix_[0][i].aW = 0.0f;
                coefMatrix_[0][i].aN = 0.0f;
                coefMatrix_[0][i].aS = 0.0f;
                coefMatrix_[0][i].bsrc = u_bc;
            }
        }
        
        // 处理 v 速度边界条件
        if (direction_ == 1)
        {
            if (bc.velocityType == WALL || bc.velocityType == INLET)
            {
                float v_bc = bc.VelocityValue[1];
                coefMatrix_[0][i].aP = 1.0f;
                coefMatrix_[0][i].aE = 0.0f;
                coefMatrix_[0][i].aW = 0.0f;
                coefMatrix_[0][i].aN = 0.0f;
                coefMatrix_[0][i].aS = 0.0f;
                coefMatrix_[0][i].bsrc = v_bc;
            }
        }
    }
    
    // === 北侧边界 (j=ncy-1) ===
    for (int i = 0; i < ncx; ++i)
    {
        const auto& bc = boundaryField_.north(i);
        
        // 处理 u 速度边界条件
        if (direction_ == 0)
        {
            if (bc.velocityType == WALL || bc.velocityType == INLET)
            {
                float u_bc = bc.VelocityValue[0];
                coefMatrix_[ncy - 1][i].aP = 1.0f;
                coefMatrix_[ncy - 1][i].aE = 0.0f;
                coefMatrix_[ncy - 1][i].aW = 0.0f;
                coefMatrix_[ncy - 1][i].aN = 0.0f;
                coefMatrix_[ncy - 1][i].aS = 0.0f;
                coefMatrix_[ncy - 1][i].bsrc = u_bc;
            }
        }
        
        // 处理 v 速度边界条件
        if (direction_ == 1)
        {
            if (bc.velocityType == WALL || bc.velocityType == INLET)
            {
                float v_bc = bc.VelocityValue[1];
                coefMatrix_[ncy - 1][i].aP = 1.0f;
                coefMatrix_[ncy - 1][i].aE = 0.0f;
                coefMatrix_[ncy - 1][i].aW = 0.0f;
                coefMatrix_[ncy - 1][i].aN = 0.0f;
                coefMatrix_[ncy - 1][i].aS = 0.0f;
                coefMatrix_[ncy - 1][i].bsrc = v_bc;
            }
        }
    }
}

void ScalarEquation::setRelaxation(float alpha)
{
    // TODO: 实现欠松弛逻辑 (Under-relaxation)
    // aP = aP / alpha;
    // bsrc = bsrc + (1 - alpha) * aP_old * scalarField_old;
}
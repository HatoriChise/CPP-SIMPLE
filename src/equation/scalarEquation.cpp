#include "scalarEquation.hpp"

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

float ScalarEquation::computeFaceMassFlux(int i, int j, int face) const
{
    auto meshSize = mesh_.getMeshSize();
    float dx = meshSize[0];
    float dy = meshSize[1];
    float rho = fluidPropertyField_(i, j).rho;

    float u_face, v_face;
    switch (face)
    {
    case 0:  // 东侧 (i+1/2, j)
        u_face = 0.5f * (vectorField_.u()(i, j) + vectorField_.u()(i + 1, j));
        v_face = 0.5f * (vectorField_.v()(i, j) + vectorField_.v()(i + 1, j));
        return rho * u_face * dy;
    case 1:  // 西侧 (i-1/2, j)
        u_face = 0.5f * (vectorField_.u()(i - 1, j) + vectorField_.u()(i, j));
        v_face = 0.5f * (vectorField_.v()(i - 1, j) + vectorField_.v()(i, j));
        return rho * u_face * dy;
    case 2:  // 北侧 (i, j+1/2)
        u_face = 0.5f * (vectorField_.u()(i, j) + vectorField_.u()(i, j + 1));
        v_face = 0.5f * (vectorField_.v()(i, j) + vectorField_.v()(i, j + 1));
        return rho * v_face * dx;
    case 3:  // 南侧 (i, j-1/2)
        u_face = 0.5f * (vectorField_.u()(i, j - 1) + vectorField_.u()(i, j));
        v_face = 0.5f * (vectorField_.v()(i, j - 1) + vectorField_.v()(i, j));
        return rho * v_face * dx;
    }
    return 0.0f;
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

void ScalarEquation::addConvectionTerm()
{
    // 遍历内部单元
    for (int j = 1; j < ncy - 1; ++j)
    {
        for (int i = 1; i < ncx - 1; ++i)
        {
            // 计算界面质量通量
            float F_e = computeFaceMassFlux(i, j, 0);  // 东 East
            float F_w = computeFaceMassFlux(i, j, 1);  // 西 West
            float F_n = computeFaceMassFlux(i, j, 2);  // 北 North
            float F_s = computeFaceMassFlux(i, j, 3);  // 南 South

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
    // TODO: 1. 获取 boundaryField_ 中的 BC 描述
    // TODO: 2. 处理 Dirichlet: 修改 aP 并将贡献转移至 bsrc (bsrc += a_nb * Phi_bc)
    // TODO: 3. 处理 Neumann: 修改 aP 且令 a_nb = 0
}

void ScalarEquation::setRelaxation(float alpha)
{
    // TODO: 实现欠松弛逻辑 (Under-relaxation)
    // aP = aP / alpha;
    // bsrc = bsrc + (1 - alpha) * aP_old * scalarField_old;
}
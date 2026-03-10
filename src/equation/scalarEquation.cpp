#include "scalarEquation.hpp"

ScalarEquation::ScalarEquation(StructuredMesh &mesh, ScalarField &scalarField,
                               VectorField &vectorField, BoundaryField &boundaryField,
                               FluidPropertyField &fluidPropertyField)
    : mesh_(mesh),
      scalarField_(scalarField),
      vectorField_(vectorField),
      boundaryField_(boundaryField),
      fluidPropertyField_(fluidPropertyField)
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
    // TODO: 1. 获取面速度 (Face Velocity) 或通过 Rhie-Chow 计算界面通量
    // TODO: 2. 根据 scheme (0: Upwind, 1: Central) 计算对流系数
    //    Upwind 示例: aE_conv = max(-F_e, 0), aW_conv = max(F_w, 0)
    // TODO: 3. 更新 aE, aW, aN, aS 及其对 aP 的贡献
}

void ScalarEquation::addSourceTerm()
{
    // TODO: 将自定义源项场 S 的值映射并累加到 bsrc: bsrc += S(i,j) * volume
}

void ScalarEquation::addPressureGradient()
{
    // TODO: 专门为动量方程设计
    // direction == 0 (U方向): bsrc += -(p_e - p_w) * dy
    // direction == 1 (V方向): bsrc += -(p_n - p_s) * dx
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
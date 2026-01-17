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
    // TODO: 1. 根据 mesh 的尺寸 (ncx, ncy) 初始化 coefMatrix_ 的大小
    // coefMatrix_.resize(boost::extents[mesh_.ncy()][mesh_.ncx()]);
    
    // TODO: 2. 执行一次 resetCoefficients() 确保初始状态为零
}

ScalarEquation::~ScalarEquation() = default;

void ScalarEquation::resetCoefficients()
{
    // TODO: 遍历 coefMatrix_，将所有单元的 aE, aW, aN, aS, aP, bsrc 重置为 0
    // 提示：这是每一轮 SIMPLE 迭代开始前必须执行的操作
}

// === 核心组装方法 (待实现) ===

void ScalarEquation::addDiffusionTerm()
{
    // TODO: 1. 获取网格几何信息 (dx, dy) 和物性 (Gamma = mu 或 k)
    // TODO: 2. 遍历内部单元 (1 to n-1)，计算界面扩散通量：
    //    aE = Gamma_e * dy / dx_Pe
    //    aW = Gamma_w * dy / dx_Pw
    //    ... 同理计算 aN, aS
    // TODO: 3. 累加到 aP: aP += (aE + aW + aN + aS)
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
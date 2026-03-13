# Rhie-Chow 插值实现计划

## 背景

当前 `computeFaceMassFlux` 使用线性插值计算界面速度，会导致棋盘格压力场问题。需要升级为 Rhie-Chow 插值。

## Rhie-Chow 公式

这里是纯数学推导：定义已知量（基于动量方程）我们在单元中心 $P$ 和 $E$ 的动量方程定义如下：

$$u_P=\hat{u}_P-d_P\left( \frac{\partial p}{\partial x} \right) _P\Longrightarrow \hat{u}_P=u_P+d_P\left( \frac{\partial p}{\partial x} \right) _P
\\
u_E=\hat{u}_E-d_E\left( \frac{\partial p}{\partial x} \right) _E\Longrightarrow \hat{u}_E=u_E+d_E\left( \frac{\partial p}{\partial x} \right) _E$$

其中 $d_P = \frac{V_P}{a_P^P}$，$d_E = \frac{V_E}{a_P^E}$。

对 $\hat{u}$ 进行线性插值我们假设界面 $e$ 处的伪速度 $\hat{u}_e$ 可以通过 $P$ 和 $E$ 处的伪速度进行加权平均得到（设权值为 $\alpha$，对于均匀网格 $\alpha=0.5$）：

$$\hat{u}_e = \alpha \hat{u}_P + (1-\alpha) \hat{u}_E$$

代入定义式（消除 $\hat{u}$ 的过程）将第 1 步中的 $\hat{u}_P$ 和 $\hat{u}_E$ 的表达式代入上式：$$\hat{u}_e = \alpha \left[ u_P + d_P \left( \frac{\partial p}{\partial x} \right)_P \right] + (1-\alpha) \left[ u_E + d_E \left( \frac{\partial p}{\partial x} \right)_E \right]$$根据界面速度的定义（界面同样满足动量平衡）：

$$u_e = \hat{u}_e - d_e \left( \frac{\partial p}{\partial x} \right)_e$$现在，把刚才展开的 $\hat{u}_e$ 代入：$$u_e = \underbrace{\left\{ \alpha u_P + (1-\alpha) u_E \right\}}_{\text{线性插值速度 } \overline{u}_e} + \underbrace{\left\{ \alpha d_P (\nabla p)_P + (1-\alpha) d_E (\nabla p)_E - d_e (\nabla p)_e \right\}}_{\text{Rhie-Chow 修正项}}$$

总结公式最终，你只需要在代码中实现：

$$u_e=\frac{1}{2}\left( u_P+u_E \right) +\frac{1}{2}d_P\left( \frac{\partial p}{\partial x} \right) _P+\frac{1}{2}d_E\left( \frac{\partial p}{\partial x} \right) _E-d_e\left( \frac{\partial p}{\partial x} \right) _e
\\
d_e=\frac{1}{2}\left( d_P+d_E \right)
\\
\left( \frac{\partial p}{\partial x} \right) _e=\left( p_P-p_E \right) $$

## 修改文件

- `src/equation/scalarEquation.hpp`
- `src/equation/scalarEquation.cpp`

## 实现方案

### 更新测试代码

## 调用顺序要求

**重要：** Rhie-Chow 插值需要动量系数 `aP`，因此调用顺序必须为：

1. `addDiffusionTerm()` - 先计算 aP
2. `addConvectionTerm(&pressure)` - 再使用 RC 插值

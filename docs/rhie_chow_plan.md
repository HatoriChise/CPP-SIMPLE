# Rhie-Chow 插值实现计划

## 背景

当前 `computeFaceMassFlux` 使用线性插值计算界面速度，会导致棋盘格压力场问题。需要升级为 Rhie-Chow 插值。

## Rhie-Chow 公式

$$u_e = \frac{1}{2}(u_P + u_E) + \frac{dy}{a_P}(p_P - p_E)$$

其中：
- $\bar{u}_e = \frac{1}{2}(u_P + u_E)$：线性插值速度
- $d_e = \frac{A_e}{a_P}$：界面动量系数倒数
- $A_e$：界面面积（东侧为 dy，北侧为 dx）
- $a_P$：动量方程中心系数

## 修改文件

- `src/equation/scalarEquation.hpp`
- `src/equation/scalarEquation.cpp`

## 实现方案

### 1. 修改方法签名（hpp）

```cpp
// 计算界面质量通量（使用Rhie-Chow插值）
// face: 0=东, 1=西, 2=北, 3=南
// pressure: 压力场指针，若为nullptr则使用线性插值
float computeFaceMassFlux(int i, int j, int face, const ScalarField* pressure) const;

void addConvectionTerm(const ScalarField* pressure = nullptr);
```

### 2. 实现 Rhie-Chow 插值（cpp）

```cpp
float ScalarEquation::computeFaceMassFlux(int i, int j, int face, const ScalarField* pressure) const
{
    auto meshSize = mesh_.getMeshSize();
    float dx = meshSize[0];
    float dy = meshSize[1];
    float rho = fluidPropertyField_(i, j).rho;

    float u_face, v_face, area;
    float dp = 0.0f;  // 压力差
    
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
    default:
        return 0.0f;
    }

    // Rhie-Chow 修正
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
```

### 3. 更新 addConvectionTerm

```cpp
void ScalarEquation::addConvectionTerm(const ScalarField* pressure)
{
    for (int j = 1; j < ncy - 1; ++j)
    {
        for (int i = 1; i < ncx - 1; ++i)
        {
            float F_e = computeFaceMassFlux(i, j, 0, pressure);
            float F_w = computeFaceMassFlux(i, j, 1, pressure);
            float F_n = computeFaceMassFlux(i, j, 2, pressure);
            float F_s = computeFaceMassFlux(i, j, 3, pressure);

            coefMatrix_[j][i].aE += std::max(-F_e, 0.0f);
            coefMatrix_[j][i].aW += std::max(F_w, 0.0f);
            coefMatrix_[j][i].aN += std::max(-F_n, 0.0f);
            coefMatrix_[j][i].aS += std::max(F_s, 0.0f);
            coefMatrix_[j][i].aP += (F_e - F_w + F_n - F_s);
        }
    }
}
```

### 4. 更新测试代码

```cpp
// 组装 u 动量方程
u_momentum.resetCoefficients();
u_momentum.addDiffusionTerm();                    // 先计算 aP
u_momentum.addConvectionTerm(&pressure);          // 传入压力场启用 RC 插值
u_momentum.addPressureGradient(pressure);
```

## 调用顺序要求

**重要：** Rhie-Chow 插值需要动量系数 `aP`，因此调用顺序必须为：

1. `addDiffusionTerm()` - 先计算 aP
2. `addConvectionTerm(&pressure)` - 再使用 RC 插值

## 后续优化

- 可以考虑将 `aP` 单独存储，支持更复杂的 RC 实现
- 可以添加松弛因子支持
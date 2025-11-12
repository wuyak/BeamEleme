# Timoshenko梁解析解代码修正

本项目使用Khasawneh 2019的Timoshenko梁解析解代码，修正了3处bug。

## 代码来源

参考代码：[Timoshenko Beam Frequencies and Modes](http://dx.doi.org/10.17632/r275tx2yp8.1)
作者：Firas A. Khasawneh, Daniel J. Segalman
论文："Exact and Numerically Stable Expressions for Euler-Bernoulli and Timoshenko Beam Modes" (2019)

我们的实现在 `src/+analytical/+bending_shear/` 目录下：
- `getEigenvalues.m` - 计算特征方程（对应原始代码的 `get_TB_eigVal.m`）
- `getModes.m` - 计算模态振型（对应原始代码的 `get_TB_eigModes.m`）

## 修正的bug

### Bug 1: CP边界条件λ=α变量使用错误

**文件**：`src/+analytical/+bending_shear/getModes.m`
**位置**：第243-244行
**边界条件**：CP (Clamped-Pinned，固定-简支)

原始代码定义了 `w2Vec = wN(alphaN)` 用于λ=α情况，但在计算系数A2、A3时错误使用了 `w1Vec`（来自λ<α区域）。

**原始代码**（错误）：
```matlab
w2Vec = wN(alphaN);
A1 = ones(1, length(lambda_EQ));
A2 = w1Vec.*sin(w1Vec)./((lambda_EQ-w1Vec.^2).*(cos(w1Vec)-1));  % 错误使用w1Vec
A3 = w1Vec./(lambda_EQ-w1Vec.^2);  % 错误使用w1Vec
```

**修正后**：
```matlab
w2Vec = wN(alphaN);
A1 = ones(1, length(lambda_EQ));
A2 = w2Vec.*sin(w2Vec)./((lambda_EQ-w2Vec.^2).*(cos(w2Vec)-1));
A3 = w2Vec./(lambda_EQ-w2Vec.^2);
```

### Bug 2: PP边界条件λ=α变量使用错误

**文件**：`src/+analytical/+bending_shear/getModes.m`
**位置**：第397-398行
**边界条件**：PP (Pinned-Pinned，简支-简支)

与Bug 1类似，在稳定表达式分支中错误使用了 `w1Vec` 而非 `w2Vec`。

**原始代码**（错误）：
```matlab
w2Vec = wN(alphaN);
eigmode_U_EQ = chebfun(@(z) sin(w1Vec.*z), z_range);  % 错误使用w1Vec
eigmode_phi_EQ = chebfun(@(z) -(lambda_EQ-w1Vec.^2)./w1Vec.*cos(w1Vec.*z), z_range);  % 错误使用w1Vec
```

**修正后**：
```matlab
w2Vec = wN(alphaN);
eigmode_U_EQ = chebfun(@(z) sin(w2Vec.*z), z_range);
eigmode_phi_EQ = chebfun(@(z) -(lambda_EQ-w2Vec.^2)./w2Vec.*cos(w2Vec.*z), z_range);
```

### Bug 3: CP边界条件缺失搜索范围

**文件**：`src/+analytical/+bending_shear/getEigenvalues.m`
**位置**：第147行
**边界条件**：CP (Clamped-Pinned，固定-简支)

原始代码在CP边界的λ>α区域定义chebfun时，遗漏了搜索范围参数 `range_GT`，导致chebfun使用默认范围 [-1, 1]，无法找到 λ>α 的特征值。

**原始代码**（错误）：
```matlab
eigVal_GT = chebfun(@(x)..., 'splitting', 'on');  % 缺少range_GT
```

**修正后**：
```matlab
eigVal_GT = chebfun(@(x)..., range_GT, 'splitting', 'on');
```

**影响**：修正前CP边界在高阶模态会缺失大量特征值。例如α=100时请求5个模态，只能找到3个（λ<α区域的全部根），缺失2个λ>α区域的模态。修正后可以正常找到所有模态。

## 验证

所有修正已通过以下方式验证：
1. 与有限元结果对比（频率误差<0.1%）
2. 与其他边界条件的代码结构对比（CC边界正确使用了w2Vec）
3. 极限情况测试（α→∞退化为Euler-Bernoulli梁）

## 引用

本项目使用的Timoshenko梁解析解算法基于以下工作：

> Firas A. Khasawneh, Daniel J. Segalman (2019).
> "Exact and Numerically Stable Expressions for Euler-Bernoulli and Timoshenko Beam Modes"
> International Journal of Mechanical Sciences, 166, 105234.
> Code: http://dx.doi.org/10.17632/r275tx2yp8.1

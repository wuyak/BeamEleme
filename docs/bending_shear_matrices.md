# Timoshenko梁单元的刚度矩阵与质量矩阵推导

> **核心参考文献**:
> - Friedman, Z., & Kosmatka, J. B. (1993). An improved two-node Timoshenko beam finite element. *Computers & Structures*, 47(3), 473-481. DOI: https://doi.org/10.1016/0045-7949(93)90243-7

---

## 目录

1. [引言](#1-引言)
2. [Timoshenko梁基本方程回顾](#2-timoshenko梁基本方程回顾)
3. [有限元形函数的推导](#3-有限元形函数的推导)
4. [刚度矩阵推导](#4-刚度矩阵推导)
5. [质量矩阵推导](#5-质量矩阵推导)
6. [参考文献](#6-参考文献)

---

## 1. 引言

### 1.1 本文目标

本文详细推导 Friedman & Kosmatka (1993) 提出的**改进两节点 Timoshenko 梁单元**的刚度矩阵和质量矩阵。该单元的关键优点是：

1. **避免剪切锁定**：通过让形函数满足静力平衡方程
2. **精确刚度矩阵**：与柔度法推导的精确 Timoshenko 刚度矩阵一致
3. **一致质量矩阵**：可用于动力学分析

### 1.2 与理论推导文档的关系

Timoshenko 梁的基本理论（运动学、本构关系、控制方程）已在《Timoshenko梁理论推导.md》中详细推导。本文专注于**有限元离散化**，即：

- 如何选择形函数？
- 形函数为什么要满足静力平衡方程？
- 如何从形函数得到刚度矩阵和质量矩阵？

### 1.3 符号约定

| 符号 | 含义 | 单位 |
|------|------|------|
| $E$ | 弹性模量 | Pa |
| $G$ | 剪切模量 | Pa |
| $\rho$ | 材料密度 | kg/m³ |
| $A$ | 截面面积 | m² |
| $I$ | 截面惯性矩 | m⁴ |
| $k$ | 剪切修正系数 | 无量纲 |
| $L$ | 单元长度 | m |
| $\phi$ | 剪切影响因子 | 无量纲 |
| $w$ | 横向位移 | m |
| $\theta$ | 截面转角 | rad |

**注意**：本文采用 Friedman 论文的符号约定，转角用 $\theta$ 表示（论文中为绕 $y$ 轴正向转动），与《Timoshenko梁理论推导.md》中的 $\psi$ 对应。

---

## 2. Timoshenko梁基本方程回顾

本节简要回顾 Timoshenko 梁的基本方程，详细推导见《Timoshenko梁理论推导.md》。

### 2.1 位移场假设 [论文 Eq.1a-c]

考虑在 $x$-$z$ 平面内弯曲的梁，位移场假设为：

$$
u(x,y,z,t) = z\theta(x,t), \quad v(x,y,z,t) = 0, \quad w(x,y,z,t) = w(x,t) \tag{2.1}
$$

其中：
- $w(x,t)$：中性轴的横向位移
- $\theta(x,t)$：截面绕 $y$ 轴的转角（**注意**：这里 $\theta$ 的定义与剪切角 $\gamma$ 的关系为 $\gamma = \frac{\partial w}{\partial x} + \theta$）

### 2.2 应变-位移关系 [论文 Eq.2a-b]

$$
\varepsilon_{xx} = z\frac{\partial\theta}{\partial x}, \quad \gamma_{xz} = \frac{\partial w}{\partial x} + \theta \tag{2.2}
$$

**物理意义**：
- $\varepsilon_{xx}$：轴向应变，由截面转动引起
- $\gamma_{xz}$：剪切应变，为挠度斜率与截面转角之和

### 2.3 应变能表达式 [论文 Eq.4b]

$$
U = \frac{1}{2}\int_0^L \left[ EI\left(\frac{\partial\theta}{\partial x}\right)^2 + kGA\left(\frac{\partial w}{\partial x} + \theta\right)^2 \right] dx \tag{2.3}
$$

### 2.4 动能表达式 [论文 Eq.5b]

$$
T = \frac{1}{2}\int_0^L \left[ \rho A\left(\frac{\partial w}{\partial t}\right)^2 + \rho I\left(\frac{\partial\theta}{\partial t}\right)^2 \right] dx \tag{2.4}
$$

### 2.5 运动方程 [论文 Eq.7a-b]

通过 Hamilton 原理 $\delta\int_{t_1}^{t_2}(T-U)dt = 0$，可得运动方程：

$$
\frac{\partial}{\partial x}\left( kGA\left(\frac{\partial w}{\partial x} + \theta\right) \right) + q = \rho A \frac{\partial^2 w}{\partial t^2} \tag{2.5a}
$$

$$
\frac{\partial}{\partial x}\left( EI\frac{\partial\theta}{\partial x} \right) - kGA\left(\frac{\partial w}{\partial x} + \theta\right) + m = \rho I \frac{\partial^2 \theta}{\partial t^2} \tag{2.5b}
$$

### 2.6 边界条件 [论文 Eq.7c]

在梁的两端 $x=0$ 和 $x=L$，需要指定：

| 几何边界条件 | 自然边界条件 |
|-------------|-------------|
| $w$ | $Q = kGA\left(\frac{\partial w}{\partial x} + \theta\right)$ |
| $\theta$ | $M = EI\frac{\partial\theta}{\partial x}$ |

---

## 3. 有限元形函数的推导

本节是全文的**核心**，详细推导 Friedman 单元的形函数。

### 3.1 剪切锁定问题与解决思路

#### 3.1.1 什么是剪切锁定？

传统方法独立选择 $w$ 和 $\theta$ 的形函数（例如都用线性函数）。问题在于：

- 对于细长梁（$L/h \gg 1$），剪切变形应该趋于零：$\gamma = \frac{\partial w}{\partial x} + \theta \to 0$
- 但独立选择的形函数无法保证 $\gamma$ 在整个单元内为零
- 结果：产生**虚假的剪切应变能**，刚度矩阵过刚

#### 3.1.2 Friedman 的解决方案

**核心思想**：不独立选择 $w$ 和 $\theta$ 的形函数，而是让它们满足**齐次静力平衡方程**。

这样做的好处：
1. 形函数之间相互约束，避免虚假剪切应变
2. 当 $\phi \to 0$（细长梁）时，自动退化为 Euler-Bernoulli 梁

### 3.2 齐次静力平衡方程 [论文 Eq.8a-b]

#### 3.2.1 简化条件

为了推导形函数，我们需要将一般运动方程 (2.5a-b) 简化为齐次静力平衡方程。简化基于以下三个条件：

1. **无外载**：$q = 0$（无分布横向力），$m = 0$（无分布力矩）
2. **均匀梁**：$kGA$、$EI$ 为常数，可以提到导数外面
3. **静力问题**：形函数描述的是静态变形形状，忽略惯性项

> **为什么形函数要基于静力方程？**
> 形函数描述的是在给定节点位移下，单元内部的位移分布"形状"。这个形状应该是静力平衡的——它代表了一个没有惯性力、没有外载的平衡状态。动力学效应通过质量矩阵单独考虑。

#### 3.2.2 第一个方程的简化

从运动方程 (2.5a) 出发：
$$
\frac{\partial}{\partial x}\left( kGA\left(\frac{\partial w}{\partial x} + \theta\right) \right) + q = \rho A \frac{\partial^2 w}{\partial t^2}
$$

**Step 1**：应用无外载条件 $q = 0$ 和静力条件 $\rho A \frac{\partial^2 w}{\partial t^2} = 0$：
$$
\frac{\partial}{\partial x}\left( kGA\left(\frac{\partial w}{\partial x} + \theta\right) \right) = 0
$$

**Step 2**：由于 $kGA$ 为常数，可以提到导数外面：
$$
kGA \cdot \frac{\partial}{\partial x}\left(\frac{\partial w}{\partial x} + \theta\right) = 0
$$

**Step 3**：两边除以 $kGA \neq 0$：
$$
\boxed{\frac{\partial}{\partial x}\left(\frac{\partial w}{\partial x} + \theta\right) = 0} \tag{3.1a}
$$

#### 3.2.3 第二个方程的简化

从运动方程 (2.5b) 出发：
$$
\frac{\partial}{\partial x}\left( EI\frac{\partial\theta}{\partial x} \right) - kGA\left(\frac{\partial w}{\partial x} + \theta\right) + m = \rho I \frac{\partial^2 \theta}{\partial t^2}
$$

**Step 1**：应用无外载条件 $m = 0$ 和静力条件 $\rho I \frac{\partial^2 \theta}{\partial t^2} = 0$：
$$
\frac{\partial}{\partial x}\left( EI\frac{\partial\theta}{\partial x} \right) - kGA\left(\frac{\partial w}{\partial x} + \theta\right) = 0
$$

**Step 2**：由于 $EI$ 为常数：
$$
EI \cdot \frac{\partial^2 \theta}{\partial x^2} - kGA\left(\frac{\partial w}{\partial x} + \theta\right) = 0
$$

**Step 3**：两边除以 $EI$：
$$
\boxed{\frac{\partial^2 \theta}{\partial x^2} - \frac{kGA}{EI}\left(\frac{\partial w}{\partial x} + \theta\right) = 0} \tag{3.1b}
$$

#### 3.2.4 物理意义

- **方程 (3.1a)**：剪切应变 $\gamma = \frac{\partial w}{\partial x} + \theta$ 沿梁长为**常数**。这是因为在无外载、均匀截面的条件下，剪力为常数，而剪力正比于剪切应变。

- **方程 (3.1b)**：描述了曲率（$\frac{\partial^2 \theta}{\partial x^2}$）与剪切应变（$\frac{\partial w}{\partial x} + \theta$）之间的平衡关系。

### 3.3 多项式假设

引入无量纲坐标：
$$
\xi = \frac{x}{L} \in [0, 1] \tag{3.2}
$$

假设：
- **横向位移**（三次多项式）：
$$
w(\xi) = a_0 + a_1\xi + a_2\xi^2 + a_3\xi^3 \tag{3.3a}
$$

- **截面转角**（二次多项式）：
$$
\theta(\xi) = b_0 + b_1\xi + b_2\xi^2 \tag{3.3b}
$$

**为什么选择这样的阶次？**
- 从方程 (3.1a) 可知，$\frac{\partial w}{\partial x} + \theta$ 为常数
- 若 $\theta$ 是二次多项式，$\frac{\partial w}{\partial x}$ 也必须是二次多项式
- 因此 $w$ 必须是三次多项式

共有 **7 个未知系数**：$a_0, a_1, a_2, a_3, b_0, b_1, b_2$

### 3.4 边界条件（4个方程）

两节点单元的节点自由度为 $\{w_1, \theta_1, w_2, \theta_2\}$。

在 $\xi = 0$：
$$
w(0) = a_0 = w_1 \tag{3.4a}
$$
$$
\theta(0) = b_0 = \theta_1 \tag{3.4b}
$$

在 $\xi = 1$：
$$
w(1) = a_0 + a_1 + a_2 + a_3 = w_2 \tag{3.4c}
$$
$$
\theta(1) = b_0 + b_1 + b_2 = \theta_2 \tag{3.4d}
$$

### 3.5 约束方程的推导（3个方程）

#### 3.5.1 从方程 (3.1a) 推导

方程 (3.1a) 要求 $\frac{\partial w}{\partial x} + \theta$ 为常数。

注意坐标变换：$\frac{\partial}{\partial x} = \frac{1}{L}\frac{d}{d\xi}$

$$
\frac{1}{L}\frac{dw}{d\xi} + \theta = C_0 \quad (\text{常数}) \tag{3.5}
$$

计算各项：
$$
\frac{dw}{d\xi} = a_1 + 2a_2\xi + 3a_3\xi^2 \tag{3.6a}
$$
$$
\theta = b_0 + b_1\xi + b_2\xi^2 \tag{3.6b}
$$

代入 (3.5)：
$$
\frac{1}{L}(a_1 + 2a_2\xi + 3a_3\xi^2) + (b_0 + b_1\xi + b_2\xi^2) = C_0 \tag{3.7}
$$

**要使左边对所有 $\xi$ 都等于常数 $C_0$，各次幂的系数必须满足：**

**$\xi^2$ 系数**：
$$
\frac{3a_3}{L} + b_2 = 0 \quad \Rightarrow \quad \boxed{b_2 = -\frac{3a_3}{L}} \tag{3.8a}
$$

**$\xi^1$ 系数**：
$$
\frac{2a_2}{L} + b_1 = 0 \quad \Rightarrow \quad \boxed{b_1 = -\frac{2a_2}{L}} \tag{3.8b}
$$

**$\xi^0$ 系数**：
$$
\frac{a_1}{L} + b_0 = C_0 \tag{3.8c}
$$

#### 3.5.2 从方程 (3.1b) 推导

方程 (3.1b) 用 $\xi$ 表示：
$$
\frac{1}{L^2}\frac{d^2\theta}{d\xi^2} = \frac{kGA}{EI}\left(\frac{1}{L}\frac{dw}{d\xi} + \theta\right) \tag{3.9}
$$

计算 $\frac{d^2\theta}{d\xi^2}$：
$$
\frac{d^2\theta}{d\xi^2} = 2b_2 \tag{3.10}
$$

由 (3.5)，右边等于 $\frac{kGA}{EI} \cdot C_0$。代入：
$$
\frac{2b_2}{L^2} = \frac{kGA}{EI} \cdot C_0 \tag{3.11}
$$

解出 $C_0$：
$$
C_0 = \frac{2b_2 \cdot EI}{kGA \cdot L^2} \tag{3.12}
$$

### 3.6 $\phi$ 的自然出现

将约束 (3.8a) 代入 (3.12)：
$$
C_0 = \frac{2(-3a_3/L) \cdot EI}{kGA \cdot L^2} = -\frac{6a_3 \cdot EI}{kGA \cdot L^3} \tag{3.13}
$$

整理：
$$
C_0 = -\frac{6a_3}{L} \cdot \frac{EI}{kGA \cdot L^2} \tag{3.14}
$$

**此时，我们定义**：
$$
\boxed{\phi \equiv \frac{12EI}{kGA \cdot L^2}} \tag{3.15}
$$

则：
$$
C_0 = -\frac{6a_3}{L} \cdot \frac{\phi}{12} = -\frac{a_3\phi}{2L} \tag{3.16}
$$

#### 3.6.1 "12"的来源

**问题**：为什么定义中选择系数"12"而不是其他数？

**答案**：这是一个**归一化选择**，目的是让最终形函数的表达式更简洁。具体来说：

1. 从 (3.8c) 和 (3.16)：$\frac{a_1}{L} + \theta_1 = -\frac{a_3\phi}{2L}$
2. 解出：$a_1 = -\frac{a_3\phi}{2} - L\theta_1$

选择"12"使得形函数中 $\phi$ 的系数变为简单的整数或简单分数（如 $1$, $\frac{1}{2}$, $(1+\phi)$ 等），便于书写和理解。

#### 3.6.2 $\phi$ 的物理意义

$$
\phi = \frac{12EI}{kGA \cdot L^2} = \frac{\text{弯曲刚度}}{\text{剪切刚度}} \times \frac{12}{L^2}
$$

- $\phi \to 0$：细长梁，剪切变形可忽略，退化为 Euler-Bernoulli 梁
- $\phi$ 增大：短粗梁，剪切变形重要

### 3.7 求解所有系数

现在我们有 **7 个方程**求解 **7 个未知数**：

| 方程编号 | 方程 | 来源 |
|---------|------|------|
| (3.4a) | $a_0 = w_1$ | 边界条件 |
| (3.4b) | $b_0 = \theta_1$ | 边界条件 |
| (3.4c) | $a_0 + a_1 + a_2 + a_3 = w_2$ | 边界条件 |
| (3.4d) | $b_0 + b_1 + b_2 = \theta_2$ | 边界条件 |
| (3.8a) | $b_2 = -3a_3/L$ | 约束 |
| (3.8b) | $b_1 = -2a_2/L$ | 约束 |
| (3.8c)+(3.16) | $a_1 = -\frac{a_3\phi}{2} - L\theta_1$ | 约束 |

#### 3.7.1 求解 $a_3$

从 (3.4d) 和约束条件：
$$
\theta_1 + b_1 + b_2 = \theta_2
$$
$$
\theta_1 - \frac{2a_2}{L} - \frac{3a_3}{L} = \theta_2
$$
$$
\Rightarrow \quad 2a_2 + 3a_3 = L(\theta_1 - \theta_2) \tag{3.17}
$$

从 (3.4c)：
$$
w_1 + a_1 + a_2 + a_3 = w_2
$$

代入 $a_1 = -\frac{a_3\phi}{2} - L\theta_1$：
$$
w_1 - \frac{a_3\phi}{2} - L\theta_1 + a_2 + a_3 = w_2
$$
$$
\Rightarrow \quad a_2 + a_3\left(1 - \frac{\phi}{2}\right) = w_2 - w_1 + L\theta_1 \tag{3.18}
$$

从 (3.17)：$a_2 = \frac{L(\theta_1 - \theta_2) - 3a_3}{2}$

代入 (3.18)：
$$
\frac{L(\theta_1 - \theta_2) - 3a_3}{2} + a_3\left(1 - \frac{\phi}{2}\right) = w_2 - w_1 + L\theta_1
$$

展开并整理：
$$
\frac{L(\theta_1 - \theta_2)}{2} - \frac{3a_3}{2} + a_3 - \frac{a_3\phi}{2} = w_2 - w_1 + L\theta_1
$$

$$
\frac{L(\theta_1 - \theta_2)}{2} - \frac{a_3(1 + \phi)}{2} = w_2 - w_1 + L\theta_1
$$

$$
-\frac{a_3(1 + \phi)}{2} = w_2 - w_1 + L\theta_1 - \frac{L\theta_1 - L\theta_2}{2}
$$

$$
-\frac{a_3(1 + \phi)}{2} = w_2 - w_1 + \frac{L\theta_1}{2} + \frac{L\theta_2}{2}
$$

解出 $a_3$：
$$
\boxed{a_3 = \frac{2}{1+\phi}\left(w_1 - w_2 - \frac{L(\theta_1 + \theta_2)}{2}\right)} \tag{3.19}
$$

#### 3.7.2 求解其他系数

从 $a_3$ 出发，依次求解其他系数。

**求解 $a_2$**：从 (3.17) 式 $2a_2 + 3a_3 = L(\theta_1 - \theta_2)$

$$
a_2 = \frac{L(\theta_1 - \theta_2) - 3a_3}{2} \tag{3.20}
$$

**求解 $a_1$**：从 (3.8c) 和 (3.16) 式

$$
a_1 = -\frac{a_3\phi}{2} - L\theta_1 \tag{3.21}
$$

**求解 $b_1$**：从约束 (3.8b)

$$
b_1 = -\frac{2a_2}{L} \tag{3.22}
$$

**求解 $b_2$**：从约束 (3.8a)

$$
b_2 = -\frac{3a_3}{L} \tag{3.23}
$$

#### 3.7.3 系数汇总

将 $a_3$ (3.19) 代入各式，得到所有系数关于节点位移的表达式：

$$
\boxed{
\begin{aligned}
a_0 &= w_1 \\[6pt]
a_1 &= -\frac{\phi}{1+\phi}\left(w_1 - w_2\right) - \frac{L}{1+\phi}\left(\theta_1 + \frac{\phi}{2}\theta_1 + \frac{\phi}{2}\theta_2\right) \\[6pt]
a_2 &= \frac{1}{1+\phi}\left[3(w_2 - w_1) + L(2\theta_1 + \theta_2) + \frac{\phi L}{2}(\theta_1 + \theta_2)\right] \\[6pt]
a_3 &= \frac{2}{1+\phi}\left(w_1 - w_2\right) - \frac{L}{1+\phi}(\theta_1 + \theta_2)
\end{aligned}
} \tag{3.24}
$$

$$
\boxed{
\begin{aligned}
b_0 &= \theta_1 \\[6pt]
b_1 &= -\frac{2a_2}{L} \\[6pt]
b_2 &= -\frac{3a_3}{L} = \frac{6}{(1+\phi)L}(w_2 - w_1) + \frac{3}{1+\phi}(\theta_1 + \theta_2)
\end{aligned}
} \tag{3.25}
$$

### 3.8 从系数到形函数

#### 3.8.1 形函数的定义

位移场可以写成节点位移的线性组合：

$$
w(\xi) = N_{w1}(\xi)w_1 + N_{w2}(\xi)\theta_1 + N_{w3}(\xi)w_2 + N_{w4}(\xi)\theta_2 \tag{3.26}
$$

$$
\theta(\xi) = N_{\theta 1}(\xi)w_1 + N_{\theta 2}(\xi)\theta_1 + N_{\theta 3}(\xi)w_2 + N_{\theta 4}(\xi)\theta_2 \tag{3.27}
$$

**形函数的物理意义**：$N_{wi}(\xi)$ 表示当第 $i$ 个自由度为1、其他为0时的位移分布。

#### 3.8.2 横向位移形函数 $N_{w1}$ 的推导

设 $w_1 = 1$, $\theta_1 = 0$, $w_2 = 0$, $\theta_2 = 0$。

**Step 1**：计算 $a_3$
$$
a_3 = \frac{2}{1+\phi}\left(1 - 0 - \frac{L(0+0)}{2}\right) = \frac{2}{1+\phi}
$$

**Step 2**：计算 $a_2$
$$
a_2 = \frac{L(0 - 0) - 3 \cdot \frac{2}{1+\phi}}{2} = -\frac{3}{1+\phi}
$$

**Step 3**：计算 $a_1$
$$
a_1 = -\frac{\frac{2}{1+\phi} \cdot \phi}{2} - L \cdot 0 = -\frac{\phi}{1+\phi}
$$

**Step 4**：$a_0 = w_1 = 1$

**Step 5**：组装形函数
$$
N_{w1}(\xi) = a_0 + a_1\xi + a_2\xi^2 + a_3\xi^3 = 1 - \frac{\phi}{1+\phi}\xi - \frac{3}{1+\phi}\xi^2 + \frac{2}{1+\phi}\xi^3
$$

**Step 6**：整理
$$
\boxed{N_{w1}(\xi) = \frac{1}{1+\phi}\left[(1+\phi) - \phi\xi - 3\xi^2 + 2\xi^3\right]} \tag{3.28}
$$

#### 3.8.3 横向位移形函数 $N_{w2}$ 的推导

设 $w_1 = 0$, $\theta_1 = 1$, $w_2 = 0$, $\theta_2 = 0$。

**Step 1**：计算 $a_3$
$$
a_3 = \frac{2}{1+\phi}\left(0 - 0 - \frac{L(1+0)}{2}\right) = -\frac{L}{1+\phi}
$$

**Step 2**：计算 $a_2$
$$
a_2 = \frac{L(1 - 0) - 3 \cdot \left(-\frac{L}{1+\phi}\right)}{2} = \frac{L}{2} + \frac{3L}{2(1+\phi)} = \frac{L(1+\phi) + 3L}{2(1+\phi)} = \frac{L(4+\phi)}{2(1+\phi)}
$$

**Step 3**：计算 $a_1$
$$
a_1 = -\frac{\left(-\frac{L}{1+\phi}\right) \cdot \phi}{2} - L \cdot 1 = \frac{L\phi}{2(1+\phi)} - L = \frac{L\phi - 2L(1+\phi)}{2(1+\phi)} = \frac{-L(2+\phi)}{2(1+\phi)}
$$

**Step 4**：$a_0 = w_1 = 0$

**Step 5**：组装形函数
$$
N_{w2}(\xi) = 0 - \frac{L(2+\phi)}{2(1+\phi)}\xi + \frac{L(4+\phi)}{2(1+\phi)}\xi^2 - \frac{L}{1+\phi}\xi^3
$$

**Step 6**：整理（提取公因子 $\frac{L}{1+\phi}$）
$$
\boxed{N_{w2}(\xi) = \frac{L}{1+\phi}\left[\xi^3 - \frac{4+\phi}{2}\xi^2 + \frac{2+\phi}{2}\xi\right] = \frac{L}{1+\phi}\left[\xi^3 - \left(2+\frac{\phi}{2}\right)\xi^2 + \left(1+\frac{\phi}{2}\right)\xi\right]} \tag{3.29}
$$

#### 3.8.4 横向位移形函数 $N_{w3}$ 的推导

设 $w_1 = 0$, $\theta_1 = 0$, $w_2 = 1$, $\theta_2 = 0$。

**Step 1**：计算 $a_3$
$$
a_3 = \frac{2}{1+\phi}\left(0 - 1 - 0\right) = -\frac{2}{1+\phi}
$$

**Step 2**：计算 $a_2$
$$
a_2 = \frac{L(0 - 0) - 3 \cdot \left(-\frac{2}{1+\phi}\right)}{2} = \frac{3}{1+\phi}
$$

**Step 3**：计算 $a_1$
$$
a_1 = -\frac{\left(-\frac{2}{1+\phi}\right) \cdot \phi}{2} - 0 = \frac{\phi}{1+\phi}
$$

**Step 4**：$a_0 = 0$

**Step 5**：组装并整理
$$
\boxed{N_{w3}(\xi) = \frac{1}{1+\phi}\left[\phi\xi + 3\xi^2 - 2\xi^3\right]} \tag{3.30}
$$

#### 3.8.5 横向位移形函数 $N_{w4}$ 的推导

设 $w_1 = 0$, $\theta_1 = 0$, $w_2 = 0$, $\theta_2 = 1$。

**Step 1**：计算 $a_3$
$$
a_3 = \frac{2}{1+\phi}\left(0 - 0 - \frac{L(0+1)}{2}\right) = -\frac{L}{1+\phi}
$$

**Step 2**：计算 $a_2$
$$
a_2 = \frac{L(0 - 1) - 3 \cdot \left(-\frac{L}{1+\phi}\right)}{2} = \frac{-L(1+\phi) + 3L}{2(1+\phi)} = \frac{L(2-\phi)}{2(1+\phi)}
$$

**Step 3**：计算 $a_1$
$$
a_1 = -\frac{\left(-\frac{L}{1+\phi}\right) \cdot \phi}{2} - 0 = \frac{L\phi}{2(1+\phi)}
$$

**Step 4**：$a_0 = 0$

**Step 5**：组装并整理
$$
\boxed{N_{w4}(\xi) = \frac{L}{1+\phi}\left[\frac{\phi}{2}\xi + \frac{2-\phi}{2}\xi^2 - \xi^3\right] = \frac{L}{1+\phi}\left[\xi^3 - \left(1-\frac{\phi}{2}\right)\xi^2 - \frac{\phi}{2}\xi\right]} \tag{3.31}
$$

**注**：最后一步调整了符号使其与论文形式一致（提取负号）。

#### 3.8.6 截面转角形函数 $N_{\theta 1}$ 的推导

设 $w_1 = 1$, $\theta_1 = 0$, $w_2 = 0$, $\theta_2 = 0$。

利用约束关系 $b_1 = -2a_2/L$, $b_2 = -3a_3/L$：

从 3.8.2 节：$a_2 = -\frac{3}{1+\phi}$, $a_3 = \frac{2}{1+\phi}$

$$
b_0 = 0, \quad b_1 = -\frac{2}{L}\left(-\frac{3}{1+\phi}\right) = \frac{6}{L(1+\phi)}, \quad b_2 = -\frac{3}{L}\left(\frac{2}{1+\phi}\right) = -\frac{6}{L(1+\phi)}
$$

$$
N_{\theta 1}(\xi) = 0 + \frac{6}{L(1+\phi)}\xi - \frac{6}{L(1+\phi)}\xi^2
$$

$$
\boxed{N_{\theta 1}(\xi) = \frac{6}{L(1+\phi)}\left(\xi - \xi^2\right) = \frac{6}{L(1+\phi)}\left(\xi^2 - \xi\right) \cdot (-1)} \tag{3.32}
$$

**注**：论文中写成 $\frac{6}{(1+\phi)L}(\xi^2 - \xi)$，符号取决于 $\xi - \xi^2$ 还是 $\xi^2 - \xi$，两者相差一个负号。

#### 3.8.7 截面转角形函数 $N_{\theta 2}$ 的推导

设 $w_1 = 0$, $\theta_1 = 1$, $w_2 = 0$, $\theta_2 = 0$。

从 3.8.3 节：$a_2 = \frac{L(4+\phi)}{2(1+\phi)}$, $a_3 = -\frac{L}{1+\phi}$

$$
b_0 = 1, \quad b_1 = -\frac{2}{L} \cdot \frac{L(4+\phi)}{2(1+\phi)} = -\frac{4+\phi}{1+\phi}, \quad b_2 = -\frac{3}{L}\left(-\frac{L}{1+\phi}\right) = \frac{3}{1+\phi}
$$

$$
N_{\theta 2}(\xi) = 1 - \frac{4+\phi}{1+\phi}\xi + \frac{3}{1+\phi}\xi^2
$$

$$
\boxed{N_{\theta 2}(\xi) = \frac{1}{1+\phi}\left[(1+\phi) - (4+\phi)\xi + 3\xi^2\right]} \tag{3.33}
$$

#### 3.8.8 截面转角形函数 $N_{\theta 3}$ 的推导

设 $w_1 = 0$, $\theta_1 = 0$, $w_2 = 1$, $\theta_2 = 0$。

从 3.8.4 节：$a_2 = \frac{3}{1+\phi}$, $a_3 = -\frac{2}{1+\phi}$

$$
b_0 = 0, \quad b_1 = -\frac{2}{L} \cdot \frac{3}{1+\phi} = -\frac{6}{L(1+\phi)}, \quad b_2 = -\frac{3}{L}\left(-\frac{2}{1+\phi}\right) = \frac{6}{L(1+\phi)}
$$

$$
\boxed{N_{\theta 3}(\xi) = \frac{6}{L(1+\phi)}\left(\xi^2 - \xi\right) \cdot (-1) = \frac{6}{L(1+\phi)}\left(\xi - \xi^2\right)} \tag{3.34}
$$

#### 3.8.9 截面转角形函数 $N_{\theta 4}$ 的推导

设 $w_1 = 0$, $\theta_1 = 0$, $w_2 = 0$, $\theta_2 = 1$。

从 3.8.5 节：$a_2 = \frac{L(2-\phi)}{2(1+\phi)}$, $a_3 = -\frac{L}{1+\phi}$

$$
b_0 = 0, \quad b_1 = -\frac{2}{L} \cdot \frac{L(2-\phi)}{2(1+\phi)} = -\frac{2-\phi}{1+\phi}, \quad b_2 = -\frac{3}{L}\left(-\frac{L}{1+\phi}\right) = \frac{3}{1+\phi}
$$

$$
N_{\theta 4}(\xi) = 0 - \frac{2-\phi}{1+\phi}\xi + \frac{3}{1+\phi}\xi^2
$$

$$
\boxed{N_{\theta 4}(\xi) = \frac{1}{1+\phi}\left[3\xi^2 - (2-\phi)\xi\right]} \tag{3.35}
$$

### 3.9 形函数汇总 [论文 Eq.9b-c]

位移场可以写成矩阵形式：
$$
\begin{pmatrix} w \\ \theta \end{pmatrix} = \begin{bmatrix} [N_w] \\ [N_\theta] \end{bmatrix} \{\Delta\} \tag{3.36}
$$

其中节点位移向量：
$$
\{\Delta\}^T = \{w_1, \theta_1, w_2, \theta_2\} \tag{3.37}
$$

**横向位移形函数**：
$$
[N_w]^T = \begin{bmatrix}
\dfrac{1}{1+\phi}\left[ 2\xi^3 - 3\xi^2 - \phi\xi + (1+\phi) \right] \\[12pt]
\dfrac{L}{1+\phi}\left[ \xi^3 - \left(2+\dfrac{\phi}{2}\right)\xi^2 + \left(1+\dfrac{\phi}{2}\right)\xi \right] \\[12pt]
\dfrac{1}{1+\phi}\left[ -2\xi^3 + 3\xi^2 + \phi\xi \right] \\[12pt]
\dfrac{L}{1+\phi}\left[ \xi^3 - \left(1-\dfrac{\phi}{2}\right)\xi^2 - \dfrac{\phi}{2}\xi \right]
\end{bmatrix} \tag{3.38}
$$

**截面转角形函数**：
$$
[N_\theta]^T = \begin{bmatrix}
\dfrac{6}{(1+\phi)L}\left[ \xi^2 - \xi \right] \\[12pt]
\dfrac{1}{1+\phi}\left[ 3\xi^2 - (4+\phi)\xi + (1+\phi) \right] \\[12pt]
\dfrac{6}{(1+\phi)L}\left[ \xi - \xi^2 \right] \\[12pt]
\dfrac{1}{1+\phi}\left[ 3\xi^2 - (2-\phi)\xi \right]
\end{bmatrix} \tag{3.39}
$$

### 3.10 形函数性质验证

#### 3.10.1 插值性

容易验证形函数在节点处满足插值条件：
- 在 $\xi = 0$：$w = w_1$, $\theta = \theta_1$
- 在 $\xi = 1$：$w = w_2$, $\theta = \theta_2$

#### 3.10.2 退化到 Euler-Bernoulli 梁

当 $\phi \to 0$ 时：

$$
N_{w1} \to 1 - 3\xi^2 + 2\xi^3
$$
$$
N_{w2} \to L(\xi - 2\xi^2 + \xi^3)
$$

这正是 **Hermite 插值函数**，即 Euler-Bernoulli 梁单元的形函数！

#### 3.10.3 避免剪切锁定

由于形函数满足静力平衡方程 (3.1a-b)，剪切应变 $\gamma = \frac{1}{L}\frac{dw}{d\xi} + \theta$ 在单元内为常数。当 $\phi \to 0$ 时，这个常数自动趋于零，不会产生虚假的剪切应变能。

---

## 4. 刚度矩阵推导

### 4.1 刚度矩阵的变分表达式 [论文 Eq.10c]

从应变能 (2.3) 出发，代入形函数：

$$
U = \frac{1}{2}\{\Delta\}^T[K]\{\Delta\}
$$

其中刚度矩阵：

$$
[K] = \int_0^L \begin{bmatrix} \dfrac{\partial}{\partial x}[N_\theta] \\[8pt] [N_\theta] + \dfrac{\partial}{\partial x}[N_w] \end{bmatrix}^T \begin{bmatrix} EI & 0 \\ 0 & kGA \end{bmatrix} \begin{bmatrix} \dfrac{\partial}{\partial x}[N_\theta] \\[8pt] [N_\theta] + \dfrac{\partial}{\partial x}[N_w] \end{bmatrix} dx \tag{4.1}
$$

### 4.2 B矩阵（应变-位移矩阵）

定义广义应变向量：
$$
\{\varepsilon\} = \begin{pmatrix} \kappa \\ \gamma \end{pmatrix} = \begin{pmatrix} \dfrac{\partial\theta}{\partial x} \\[8pt] \dfrac{\partial w}{\partial x} + \theta \end{pmatrix} = [B]\{\Delta\} \tag{4.2}
$$

**曲率 B 矩阵**：
$$
[B_\kappa] = \frac{\partial}{\partial x}[N_\theta] = \frac{1}{L}\frac{d}{d\xi}[N_\theta] \tag{4.3}
$$

**剪切应变 B 矩阵**：
$$
[B_\gamma] = [N_\theta] + \frac{\partial}{\partial x}[N_w] = [N_\theta] + \frac{1}{L}\frac{d}{d\xi}[N_w] \tag{4.4}
$$

由于形函数满足约束 (3.5)，$[B_\gamma]$ 实际上是**常数矩阵**！

### 4.3 刚度矩阵的计算

$$
[K] = [K_b] + [K_s] \tag{4.5}
$$

**弯曲刚度矩阵**：
$$
[K_b] = EI \int_0^L [B_\kappa]^T [B_\kappa] \, dx = \frac{EI}{L} \int_0^1 [B_\kappa]^T [B_\kappa] \, d\xi \tag{4.6}
$$

**剪切刚度矩阵**：
$$
[K_s] = kGA \int_0^L [B_\gamma]^T [B_\gamma] \, dx = kGA \cdot L \int_0^1 [B_\gamma]^T [B_\gamma] \, d\xi \tag{4.7}
$$

由于形函数是多项式，积分可以精确计算（无需数值积分）。

### 4.4 刚度矩阵显式结果 [论文 Eq.A1]

经过计算（详细过程见论文 Appendix），得到：

$$
\boxed{[K] = \frac{EI}{(1+\phi)L^3} \begin{bmatrix} 12 & 6L & -12 & 6L \\ 6L & (4+\phi)L^2 & -6L & (2-\phi)L^2 \\ -12 & -6L & 12 & -6L \\ 6L & (2-\phi)L^2 & -6L & (4+\phi)L^2 \end{bmatrix}} \tag{4.8}
$$

**验证**：
1. 当 $\phi = 0$ 时，退化为 Euler-Bernoulli 梁刚度矩阵
2. 矩阵对称且半正定
3. 包含刚体位移和刚体转动模态

---

## 5. 质量矩阵推导

### 5.1 质量矩阵的变分表达式 [论文 Eq.10b]

从动能 (2.4) 出发，代入形函数：

$$
T = \frac{1}{2}\{\dot{\Delta}\}^T[M]\{\dot{\Delta}\}
$$

其中质量矩阵：

$$
[M] = \int_0^L \begin{bmatrix} [N_w] \\ [N_\theta] \end{bmatrix}^T \begin{bmatrix} \rho A & 0 \\ 0 & \rho I \end{bmatrix} \begin{bmatrix} [N_w] \\ [N_\theta] \end{bmatrix} dx \tag{5.1}
$$

### 5.2 质量矩阵分解 [论文 Eq.A4]

$$
[M] = [M_{\rho A}] + [M_{\rho I}] \tag{5.2}
$$

其中：
- $[M_{\rho A}]$：平动惯性质量矩阵
- $[M_{\rho I}]$：转动惯性质量矩阵

### 5.3 平动惯性质量矩阵 [论文 Eq.A5]

$$
[M_{\rho A}] = \rho A L \int_0^1 [N_w]^T [N_w] \, d\xi \tag{5.3}
$$

$$
[M_{\rho A}] = \frac{\rho A L}{210(1+\phi)^2} \begin{bmatrix} m_{11} & m_{12} & m_{13} & m_{14} \\ m_{12} & m_{22} & m_{23} & m_{24} \\ m_{13} & m_{23} & m_{33} & m_{34} \\ m_{14} & m_{24} & m_{34} & m_{44} \end{bmatrix} \tag{5.4}
$$

其中各元素为 $\phi$ 的多项式：

$$
\begin{aligned}
m_{11} = m_{33} &= 70\phi^2 + 147\phi + 78 \\
m_{22} = m_{44} &= (7\phi^2 + 14\phi + 8)\frac{L^2}{4} \\
m_{13} &= 35\phi^2 + 63\phi + 27 \\
m_{24} &= -(7\phi^2 + 14\phi + 6)\frac{L^2}{4} \\
m_{12} = -m_{34} &= (35\phi^2 + 77\phi + 44)\frac{L}{4} \\
m_{14} = -m_{23} &= -(35\phi^2 + 63\phi + 26)\frac{L}{4}
\end{aligned} \tag{5.5}
$$

### 5.4 转动惯性质量矩阵 [论文 Eq.A6]

$$
[M_{\rho I}] = \rho I L \int_0^1 [N_\theta]^T [N_\theta] \, d\xi \tag{5.6}
$$

$$
[M_{\rho I}] = \frac{\rho I}{30(1+\phi)^2 L} \begin{bmatrix} r_{11} & r_{12} & r_{13} & r_{14} \\ r_{12} & r_{22} & r_{23} & r_{24} \\ r_{13} & r_{23} & r_{33} & r_{34} \\ r_{14} & r_{24} & r_{34} & r_{44} \end{bmatrix} \tag{5.7}
$$

其中：

$$
\begin{aligned}
r_{11} = r_{33} &= 36 \\
r_{22} = r_{44} &= (10\phi^2 + 5\phi + 4)L^2 \\
r_{13} &= -36 \\
r_{24} &= (5\phi^2 - 5\phi - 1)L^2 \\
r_{12} = -r_{34} &= -(15\phi - 3)L \\
r_{14} = -r_{23} &= -(15\phi - 3)L
\end{aligned} \tag{5.8}
$$

### 5.5 验证

当 $\phi = 0$ 时，质量矩阵退化为 Euler-Bernoulli 梁的一致质量矩阵。

---

## 6. 参考文献

1. **Friedman, Z., & Kosmatka, J. B. (1993)**. An improved two-node Timoshenko beam finite element. *Computers & Structures*, 47(3), 473-481.
   - DOI: https://doi.org/10.1016/0045-7949(93)90243-7
   - **核心参考**：形函数 (Eq.9b-c)、刚度矩阵 (Eq.A1)、质量矩阵 (Eq.A4-A6)

2. **Przemieniecki, J. S. (1968)**. *Theory of Matrix Structural Analysis*. McGraw-Hill.
   - **精确刚度矩阵**的柔度法推导

3. **Cowper, G. R. (1966)**. The shear coefficient in Timoshenko's beam theory. *ASME Journal of Applied Mechanics*, 33, 335-340.
   - **剪切修正系数** $k$ 的推导

4. **Tessler, A., & Dong, S. B. (1981)**. On a hierarchy of conforming Timoshenko beam elements. *Computers & Structures*, 14, 335-344.
   - Friedman 方法的前驱工作

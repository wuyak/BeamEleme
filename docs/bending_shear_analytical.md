# 弯曲-剪切模态的解析解推导

> **参考文献**:
> - Khasawneh, F. A., & Segalman, D. J. (2019). *Exact and numerically stable expressions for Euler-Bernoulli and Timoshenko beam modes*. Applied Acoustics, 141, 371-373.
>   - 传统特征方程：Table 7（shaded cells）
>   - 稳定化特征方程：Table 7, 8
>   - 模态振型：Eq. 16-18
>   - Roller边界条件：Appendix C

---

## 目录

1. [Timoshenko梁控制方程](#1-timoshenko梁控制方程)
2. [无量纲化处理](#2-无量纲化处理)
3. [分离变量法求解](#3-分离变量法求解)
4. [特征根的三区域分解](#4-特征根的三区域分解)
5. [模态振型通解](#5-模态振型通解)
6. [边界条件与特征方程](#6-边界条件与特征方程)
7. [数值稳定化](#7-数值稳定化)

---

## 1. Timoshenko梁控制方程

> **推导来源**：详见 [`docs/timoshenko_beam_theory.md`](docs/timoshenko_beam_theory.md) 第3-4章。

### 1.1 符号定义

| 符号 | 含义 | 单位 |
|-----|------|-----|
| $u(x,t)$ | 横向位移 | m |
| $\phi(x,t)$ | 截面转角（由弯曲产生） | rad |
| $\rho$ | 材料密度 | kg/m³ |
| $E$ | 弹性模量 | Pa |
| $G$ | 剪切模量 | Pa |
| $A$ | 截面面积 | m² |
| $I$ | 截面惯性矩 | m⁴ |
| $\kappa$ | 剪切修正系数 | - |
| $L$ | 梁长度 | m |

### 1.2 本构关系

**剪力**（由剪切应变产生）：
$$
V = AG\kappa\left(\frac{\partial u}{\partial x} - \phi\right) \tag{1.1}
$$

**弯矩**（由弯曲曲率产生）：
$$
M = EI\frac{\partial \phi}{\partial x} \tag{1.2}
$$

### 1.3 控制方程

Timoshenko梁的耦合控制方程：

$$
\boxed{\rho A\frac{\partial^2 u}{\partial t^2} = \frac{\partial}{\partial x}\left[AG\kappa\left(\frac{\partial u}{\partial x} - \phi\right)\right]} \tag{1.3}
$$

$$
\boxed{\rho I\frac{\partial^2 \phi}{\partial t^2} = AG\kappa\left(\frac{\partial u}{\partial x} - \phi\right) + \frac{\partial}{\partial x}\left(EI\frac{\partial \phi}{\partial x}\right)} \tag{1.4}
$$

这是**两个二阶耦合偏微分方程**，有两个自由度：横向位移 $u$ 和截面转角 $\phi$。

---

## 2. 无量纲化处理

### 2.1 无量纲变量

定义无量纲坐标：
$$
\zeta = \frac{x}{L}, \quad \tau = \frac{t}{T_T} \tag{2.1}
$$

定义无量纲位移和转角：
$$
u^*(\zeta,\tau) = \frac{u(\zeta L, \tau T_T)}{L}, \quad \phi^*(\zeta,\tau) = \phi(\zeta L, \tau T_T) \tag{2.2}
$$

其中时间尺度选取为：
$$
\boxed{T_T = L\sqrt{\frac{\rho}{G\kappa}}} \tag{2.3}
$$

这个时间尺度与**剪切波速**相关：$c_s = \sqrt{G\kappa/\rho}$，故 $T_T = L/c_s$。

### 2.2 无量纲参数

定义三个关键无量纲参数：

$$
\boxed{\alpha = \frac{AL^2}{I}} \tag{2.4}
$$

$\alpha$ 反映梁的**细长程度**。对于矩形截面 $b \times h$：
$$
\alpha = \frac{bh \cdot L^2}{bh^3/12} = \frac{12L^2}{h^2}
$$

细长梁 $L \gg h$ 时 $\alpha \gg 1$；短粗梁时 $\alpha$ 较小。

$$
\boxed{\beta = \frac{AG\kappa L^2}{EI}} \tag{2.5}
$$

$\beta$ 反映**剪切刚度**与**弯曲刚度**的比值。

$$
\boxed{\gamma = \frac{\beta}{\alpha} = \frac{G\kappa}{E}} \tag{2.6}
$$

$\gamma$ 是**纯材料参数**，与几何形状无关。对于各向同性材料：
$$
\gamma = \frac{G\kappa}{E} = \frac{\kappa}{2(1+\nu)}
$$

其中 $\nu$ 为泊松比。对于钢材（$\nu \approx 0.3$，$\kappa \approx 5/6$），$\gamma \approx 0.32$。

### 2.3 无量纲控制方程

将无量纲变量代入控制方程 (1.5)、(1.6)。

**推导过程**：

对于方程 (1.5)，左边：
$$
\rho A\frac{\partial^2 u}{\partial t^2} = \rho A \cdot \frac{L}{T_T^2} \cdot \frac{\partial^2 u^*}{\partial \tau^2} = \rho A \cdot \frac{L \cdot G\kappa}{\rho L^2} \cdot \frac{\partial^2 u^*}{\partial \tau^2} = \frac{AG\kappa}{L}\frac{\partial^2 u^*}{\partial \tau^2}
$$

右边：
$$
\frac{\partial}{\partial x}\left[AG\kappa\left(\frac{\partial u}{\partial x} - \phi\right)\right] = AG\kappa \cdot \frac{1}{L^2} \cdot \left(\frac{\partial^2 u^*}{\partial \zeta^2} - \frac{\partial \phi^*}{\partial \zeta}\right) = \frac{AG\kappa}{L^2}\left(\frac{\partial^2 u^*}{\partial \zeta^2} - \frac{\partial \phi^*}{\partial \zeta}\right)
$$

两边同除以 $AG\kappa/L$：
$$
\frac{\partial^2 u^*}{\partial \tau^2} = \frac{\partial^2 u^*}{\partial \zeta^2} - \frac{\partial \phi^*}{\partial \zeta} \tag{2.7}
$$

类似地，对方程 (1.6) 进行无量纲化：
$$
\frac{\partial^2 \phi^*}{\partial \tau^2} = \frac{1}{\gamma}\frac{\partial^2 \phi^*}{\partial \zeta^2} + \alpha\frac{\partial u^*}{\partial \zeta} - \alpha\phi^* \tag{2.8}
$$

**无量纲控制方程**：

$$
\boxed{\frac{\partial^2 u^*}{\partial \tau^2} = \frac{\partial^2 u^*}{\partial \zeta^2} - \frac{\partial \phi^*}{\partial \zeta}} \tag{2.9a}
$$

$$
\boxed{\frac{\partial^2 \phi^*}{\partial \tau^2} = \frac{1}{\gamma}\frac{\partial^2 \phi^*}{\partial \zeta^2} + \alpha\frac{\partial u^*}{\partial \zeta} - \alpha\phi^*} \tag{2.9b}
$$

---

## 3. 分离变量法求解

### 3.1 分离变量假设

假设解的形式为：
$$
u^*(\zeta,\tau) = U^*(\zeta)Y(\tau), \quad \phi^*(\zeta,\tau) = \Phi^*(\zeta)Y(\tau) \tag{3.1}
$$

代入无量纲控制方程 (2.9)：
$$
U^* \ddot{Y} = {U^*}'' Y - {\Phi^*}' Y
$$
$$
\Phi^* \ddot{Y} = \frac{1}{\gamma}{\Phi^*}'' Y + \alpha{U^*}' Y - \alpha\Phi^* Y
$$

分离变量：
$$
\frac{\ddot{Y}}{Y} = \frac{{U^*}'' - {\Phi^*}'}{U^*} = -\lambda \tag{3.2}
$$

取为负常数 $-\lambda$（$\lambda > 0$ 保证振动解）。

### 3.2 时间方程

$$
\ddot{Y} + \lambda Y = 0 \tag{3.3}
$$

通解为：
$$
Y(\tau) = C\cos(\sqrt{\lambda}\tau) + D\sin(\sqrt{\lambda}\tau) \tag{3.4}
$$

**无量纲角频率**为 $\sqrt{\lambda}$。

### 3.3 物理频率与特征值的关系

由 $\tau = t/T_T$ 和 $\omega = \sqrt{\lambda}/T_T$：
$$
\boxed{\hat{\omega} = \frac{\sqrt{\lambda}}{T_T} = \frac{\sqrt{\lambda}}{L}\sqrt{\frac{G\kappa}{\rho}}} \tag{3.5}
$$

物理频率（Hz）为：
$$
\boxed{f = \frac{\hat{\omega}}{2\pi} = \frac{\sqrt{\lambda}}{2\pi L}\sqrt{\frac{G\kappa}{\rho}}} \tag{3.6}
$$

### 3.4 空间方程（特征值问题）

将分离变量形式代入空间方程：
$$
-{U^*}'' + {\Phi^*}' = \lambda U^* \tag{3.7a}
$$
$$
-\frac{1}{\gamma}{\Phi^*}'' - \alpha{U^*}' + \alpha\Phi^* = \lambda\Phi^* \tag{3.7b}
$$

这是关于 $\begin{bmatrix} U^* \\ \Phi^* \end{bmatrix}$ 的**特征值问题**。

---

## 4. 特征根的三区域分解

### 4.1 指数形式假设

假设特征函数形式为：
$$
\begin{bmatrix} U^*(\zeta) \\ \Phi^*(\zeta) \end{bmatrix} = e^{m\zeta}\begin{bmatrix} w_1 \\ w_2 \end{bmatrix} \tag{4.1}
$$

代入特征值问题 (3.7)：
$$
\begin{bmatrix} -m^2 - \lambda & m \\ -\alpha m & -\frac{1}{\gamma}m^2 + (\alpha - \lambda) \end{bmatrix} \begin{bmatrix} w_1 \\ w_2 \end{bmatrix} = \begin{bmatrix} 0 \\ 0 \end{bmatrix} \tag{4.2}
$$

### 4.2 特征方程推导

非平凡解要求系数矩阵行列式为零：
$$
\det\begin{bmatrix} -m^2 - \lambda & m \\ -\alpha m & -\frac{1}{\gamma}m^2 + (\alpha - \lambda) \end{bmatrix} = 0
$$

**Step 1**：展开行列式
$$
(-m^2 - \lambda)\left(-\frac{1}{\gamma}m^2 + \alpha - \lambda\right) - m \cdot (-\alpha m) = 0
$$

**Step 2**：展开乘积
$$
\frac{1}{\gamma}m^4 - (\alpha - \lambda)m^2 + \frac{\lambda}{\gamma}m^2 - \lambda(\alpha - \lambda) + \alpha m^2 = 0
$$

**Step 3**：整理
$$
\frac{1}{\gamma}m^4 + m^2\left[-\alpha + \lambda + \frac{\lambda}{\gamma} + \alpha\right] - \lambda(\alpha - \lambda) = 0
$$
$$
\frac{1}{\gamma}m^4 + m^2 \cdot \lambda\left(1 + \frac{1}{\gamma}\right) + \lambda(\lambda - \alpha) = 0
$$

**Step 4**：乘以 $\gamma$
$$
\boxed{m^4 + \lambda(1+\gamma)m^2 + \gamma\lambda(\lambda - \alpha) = 0} \tag{4.3}
$$

### 4.3 求解 $m^2$

方程 (4.3) 是关于 $m^2$ 的二次方程。

设 $X = m^2$：
$$
X^2 + \lambda(1+\gamma)X + \gamma\lambda(\lambda - \alpha) = 0
$$

由求根公式：
$$
X = \frac{-\lambda(1+\gamma) \pm \sqrt{\lambda^2(1+\gamma)^2 - 4\gamma\lambda(\lambda - \alpha)}}{2}
$$

**Step 1**：简化判别式
$$
\Delta' = \lambda^2(1+\gamma)^2 - 4\gamma\lambda(\lambda - \alpha)
$$
$$
= \lambda^2(1+\gamma)^2 - 4\gamma\lambda^2 + 4\gamma\lambda\alpha
$$
$$
= \lambda^2\left[(1+\gamma)^2 - 4\gamma + \frac{4\gamma\alpha}{\lambda}\right]
$$
$$
= \lambda^2(1+\gamma)^2\left[1 - \frac{4\gamma}{(1+\gamma)^2} + \frac{4\gamma\alpha}{\lambda(1+\gamma)^2}\right]
$$
$$
= \lambda^2(1+\gamma)^2\left[1 - \frac{4\gamma}{(1+\gamma)^2}\left(1 - \frac{\alpha}{\lambda}\right)\right]
$$

**Step 2**：定义判别式参数
$$
\boxed{\Delta = 1 - \frac{4\gamma}{(1+\gamma)^2}\left(1 - \frac{\alpha}{\lambda}\right)} \tag{4.4}
$$

则：
$$
\sqrt{\Delta'} = \lambda(1+\gamma)\sqrt{\Delta}
$$

**Step 3**：$m^2$ 的两个解
$$
\boxed{m^2 = -\frac{1}{2}\lambda(1+\gamma)(1 \pm \sqrt{\Delta})} \tag{4.5}
$$

### 4.4 判别式 $\Delta$ 的分析

由 (4.4)，$\Delta$ 的性质取决于 $\lambda$ 与 $\alpha$ 的关系：

**Case 1**: 当 $\lambda < \alpha$ 时
$$
1 - \frac{\alpha}{\lambda} < 0 \quad \Rightarrow \quad \Delta > 1
$$

**Case 2**: 当 $\lambda = \alpha$ 时
$$
1 - \frac{\alpha}{\lambda} = 0 \quad \Rightarrow \quad \Delta = 1
$$

**Case 3**: 当 $\lambda > \alpha$ 时
$$
0 < 1 - \frac{\alpha}{\lambda} < 1 \quad \Rightarrow \quad 0 < \Delta < 1
$$

### 4.5 定义辅助变量

根据 $\Delta$ 的取值范围，定义：

**振荡率 $\omega$**（对所有情况）：
$$
\boxed{\omega^2 = \frac{1}{2}\lambda(1+\gamma)(\sqrt{\Delta} + 1)} \tag{4.6}
$$

**衰减率 $\mu$**（当 $\lambda < \alpha$，即 $\sqrt{\Delta} > 1$）：
$$
\boxed{\mu^2 = \frac{1}{2}\lambda(1+\gamma)(\sqrt{\Delta} - 1)} \tag{4.7}
$$

**高频率 $\theta$**（当 $\lambda > \alpha$，即 $\sqrt{\Delta} < 1$）：
$$
\boxed{\theta^2 = \frac{1}{2}\lambda(1+\gamma)(1 - \sqrt{\Delta})} \tag{4.8}
$$

### 4.6 特征根总结

| 区域 | 条件 | $\Delta$ | 特征根 $m$ |
|-----|------|---------|-----------|
| Region I | $\lambda < \alpha$ | $> 1$ | $\pm\mu$（实根），$\pm\omega i$（虚根） |
| Region II | $\lambda = \alpha$ | $= 1$ | $0$（二重根），$\pm\omega i$（虚根） |
| Region III | $\lambda > \alpha$ | $< 1$ | $\pm\theta i$（虚根），$\pm\omega i$（虚根） |

---

## 5. 模态振型通解

### 5.1 振型系数关系

由矩阵方程 (4.2) 的第一行：
$$
(-m^2 - \lambda)w_1 + m \cdot w_2 = 0
$$
$$
\frac{w_2}{w_1} = \frac{m^2 + \lambda}{m} = \frac{\lambda + m^2}{m} \tag{5.1}
$$

### 5.2 Region I：$\lambda < \alpha$

特征根为 $m = \pm\mu, \pm\omega i$。

**推导 $\mu$ 对应的振型分量**：

设 $m = \mu$（实根），由 (5.1)：
$$
\frac{w_2}{w_1} = \frac{\lambda + \mu^2}{\mu}
$$

对应的振型分量为：
$$
e^{\mu\zeta}\begin{bmatrix} 1 \\ \frac{\lambda + \mu^2}{\mu} \end{bmatrix}
$$

设 $m = -\mu$，由 (5.1)：
$$
\frac{w_2}{w_1} = \frac{\lambda + \mu^2}{-\mu}
$$

对应的振型分量为：
$$
e^{-\mu\zeta}\begin{bmatrix} 1 \\ -\frac{\lambda + \mu^2}{\mu} \end{bmatrix}
$$

**Step 1**：组合 $e^{\mu\zeta}$ 和 $e^{-\mu\zeta}$

设 $A_1, A_2$ 为待定系数，组合为：
$$
A_1 e^{\mu\zeta}\begin{bmatrix} 1 \\ \frac{\lambda + \mu^2}{\mu} \end{bmatrix} + A_2 e^{-\mu\zeta}\begin{bmatrix} 1 \\ -\frac{\lambda + \mu^2}{\mu} \end{bmatrix}
$$

**Step 2**：转换为双曲函数形式

令 $A_1 = (G_1 + G_2)/2$，$A_2 = (G_1 - G_2)/2$，利用：
$$
\sinh(\mu\zeta) = \frac{e^{\mu\zeta} - e^{-\mu\zeta}}{2}, \quad \cosh(\mu\zeta) = \frac{e^{\mu\zeta} + e^{-\mu\zeta}}{2}
$$

得到双曲函数形式：
$$
G_1\begin{bmatrix} \sinh(\mu\zeta) \\ \frac{\lambda + \mu^2}{\mu}\cosh(\mu\zeta) \end{bmatrix} + G_2\begin{bmatrix} \cosh(\mu\zeta) \\ \frac{\lambda + \mu^2}{\mu}\sinh(\mu\zeta) \end{bmatrix}
$$

**推导 $\omega i$ 对应的振型分量**：

设 $m = \omega i$（虚根），由 (5.1)：
$$
\frac{w_2}{w_1} = \frac{\lambda + (\omega i)^2}{\omega i} = \frac{\lambda - \omega^2}{\omega i} = -\frac{\lambda - \omega^2}{\omega}i
$$

利用 $e^{i\omega\zeta} = \cos(\omega\zeta) + i\sin(\omega\zeta)$，实部和虚部分离后得到三角函数形式：
$$
G_3\begin{bmatrix} \sin(\omega\zeta) \\ -\frac{\lambda - \omega^2}{\omega}\cos(\omega\zeta) \end{bmatrix} + G_4\begin{bmatrix} \cos(\omega\zeta) \\ \frac{\lambda - \omega^2}{\omega}\sin(\omega\zeta) \end{bmatrix}
$$

**Region I 的通解**：

$$
\boxed{\begin{bmatrix} U^*(\zeta) \\ \Phi^*(\zeta) \end{bmatrix} = A_1\begin{bmatrix} \sinh\mu\zeta \\ \frac{\lambda+\mu^2}{\mu}\cosh\mu\zeta \end{bmatrix} + A_2\begin{bmatrix} \cosh\mu\zeta \\ \frac{\lambda+\mu^2}{\mu}\sinh\mu\zeta \end{bmatrix} + A_3\begin{bmatrix} \sin\omega\zeta \\ -\frac{\lambda-\omega^2}{\omega}\cos\omega\zeta \end{bmatrix} + A_4\begin{bmatrix} \cos\omega\zeta \\ \frac{\lambda-\omega^2}{\omega}\sin\omega\zeta \end{bmatrix}} \tag{5.2}
$$

### 5.3 Region II：$\lambda = \alpha$

当 $\lambda = \alpha$ 时，判别式 $\Delta = 1$。代入 (4.5)：
$$
m^2 = -\frac{1}{2}\lambda(1+\gamma)(1 \pm 1)
$$

- 取负号：$m^2 = 0 \Rightarrow m = 0$（**二重根**）
- 取正号：$m^2 = -\lambda(1+\gamma) < 0 \Rightarrow m = \pm\omega i$

其中 $\omega^2 = \lambda(1+\gamma)$。

#### 5.3.1 二重根 $m = 0$ 的特殊处理

对于二重根，需要引入**广义特征向量**。

**Step 1**：求零空间的基向量

当 $m = 0$ 且 $\lambda = \alpha$ 时，矩阵方程 (4.2) 退化为：
$$
\begin{bmatrix} -\alpha & 0 \\ 0 & 0 \end{bmatrix}\begin{bmatrix} w_1 \\ w_2 \end{bmatrix} = \begin{bmatrix} 0 \\ 0 \end{bmatrix}
$$

零空间的基向量为：
$$
\mathbf{v}_1 = \begin{bmatrix} 0 \\ 1 \end{bmatrix}
$$

对应的振型分量为 $e^{0\cdot\zeta}\mathbf{v}_1 = \begin{bmatrix} 0 \\ 1 \end{bmatrix}$（常数转角）。

**Step 2**：求广义特征向量

由于是二重根，需要求解广义特征向量 $\mathbf{v}_2$：
$$
\mathbf{A}\mathbf{v}_2 = \mathbf{v}_1, \quad \text{其中 } \mathbf{A} = \begin{bmatrix} -\alpha & 0 \\ 0 & 0 \end{bmatrix}
$$

即：
$$
\begin{bmatrix} -\alpha & 0 \\ 0 & 0 \end{bmatrix}\begin{bmatrix} v_{21} \\ v_{22} \end{bmatrix} = \begin{bmatrix} 0 \\ 1 \end{bmatrix}
$$

第一行给出 $-\alpha v_{21} = 0 \Rightarrow v_{21} = 0$，第二行 $0 = 1$ 矛盾。

**修正**：采用标准的广义特征向量方法，对应的振型分量应为 $\zeta e^{0\cdot\zeta}\mathbf{v}_1 + e^{0\cdot\zeta}\mathbf{v}_2$。

由微分方程的通解结构，二重根 $m = 0$ 对应的两个线性无关解为：
$$
\begin{bmatrix} U^* \\ \Phi^* \end{bmatrix} = \begin{bmatrix} 0 \\ 1 \end{bmatrix}, \quad \begin{bmatrix} 1 \\ \alpha\zeta \end{bmatrix}
$$

> **验证**：代入原微分方程 (3.7)，当 $\lambda = \alpha$ 时：
> - 对于 $\begin{bmatrix} 0 \\ 1 \end{bmatrix}$：${U^*}'' = 0$，${\Phi^*}' = 0$，$\lambda U^* = 0$，方程成立
> - 对于 $\begin{bmatrix} 1 \\ \alpha\zeta \end{bmatrix}$：${U^*}'' = 0$，${\Phi^*}' = \alpha$，$-0 + \alpha = \alpha \cdot 1$，方程成立

#### 5.3.2 虚根 $\pm\omega i$ 的振型分量

与 Region I 中的推导相同，虚根 $m = \pm\omega i$ 对应的振型分量为：
$$
\begin{bmatrix} \sin\omega\zeta \\ -\frac{\lambda-\omega^2}{\omega}\cos\omega\zeta \end{bmatrix}, \quad \begin{bmatrix} \cos\omega\zeta \\ \frac{\lambda-\omega^2}{\omega}\sin\omega\zeta \end{bmatrix}
$$

#### 5.3.3 Region II 的通解

$$
\boxed{\begin{bmatrix} U^*(\zeta) \\ \Phi^*(\zeta) \end{bmatrix} = A_1\begin{bmatrix} 0 \\ 1 \end{bmatrix} + A_2\begin{bmatrix} 1 \\ \alpha\zeta \end{bmatrix} + A_3\begin{bmatrix} \sin\omega\zeta \\ -\frac{\lambda-\omega^2}{\omega}\cos\omega\zeta \end{bmatrix} + A_4\begin{bmatrix} \cos\omega\zeta \\ \frac{\lambda-\omega^2}{\omega}\sin\omega\zeta \end{bmatrix}} \tag{5.3}
$$

### 5.4 Region III：$\lambda > \alpha$

当 $\lambda > \alpha$ 时，$0 < \Delta < 1$，特征根全部为纯虚数。

#### 5.4.1 特征根分析

由 (4.4) 式：
$$
m^2 = -\frac{1}{2}\lambda(1+\gamma)(1 \pm \sqrt{\Delta})
$$

由于 $0 < \Delta < 1$，有 $0 < \sqrt{\Delta} < 1$，因此：
- $1 + \sqrt{\Delta} > 0$ → $m^2 = -\frac{1}{2}\lambda(1+\gamma)(1 + \sqrt{\Delta}) < 0$
- $1 - \sqrt{\Delta} > 0$ → $m^2 = -\frac{1}{2}\lambda(1+\gamma)(1 - \sqrt{\Delta}) < 0$

两个 $m^2$ 都为负，对应四个纯虚根 $m = \pm\omega i$ 和 $m = \pm\theta i$，其中：

$$
\omega^2 = \frac{1}{2}\lambda(1+\gamma)(1 + \sqrt{\Delta}), \quad \theta^2 = \frac{1}{2}\lambda(1+\gamma)(1 - \sqrt{\Delta})
$$

注意 $\omega > \theta > 0$（因为 $1 + \sqrt{\Delta} > 1 - \sqrt{\Delta}$）。

#### 5.4.2 $\pm\theta i$ 根的振型分量推导

设 $m = \theta i$（其中 $\theta > 0$），由 (5.1) 式：
$$
\frac{w_2}{w_1} = \frac{\lambda + m^2}{m} = \frac{\lambda + (\theta i)^2}{\theta i} = \frac{\lambda - \theta^2}{\theta i} = -\frac{\lambda - \theta^2}{\theta}i
$$

取 $w_1 = 1$，则 $w_2 = -\frac{\lambda - \theta^2}{\theta}i$。

复数形式的解为：
$$
\begin{bmatrix} U^* \\ \Phi^* \end{bmatrix} = \begin{bmatrix} 1 \\ -\frac{\lambda - \theta^2}{\theta}i \end{bmatrix} e^{i\theta\zeta}
$$

利用欧拉公式 $e^{i\theta\zeta} = \cos(\theta\zeta) + i\sin(\theta\zeta)$：
$$
\begin{bmatrix} U^* \\ \Phi^* \end{bmatrix} = \begin{bmatrix} 1 \\ -\frac{\lambda - \theta^2}{\theta}i \end{bmatrix} [\cos(\theta\zeta) + i\sin(\theta\zeta)]
$$

展开：
$$
= \begin{bmatrix} \cos(\theta\zeta) + i\sin(\theta\zeta) \\ -\frac{\lambda - \theta^2}{\theta}i\cos(\theta\zeta) + \frac{\lambda - \theta^2}{\theta}\sin(\theta\zeta) \end{bmatrix}
$$

分离实部和虚部：
- **实部**：$\begin{bmatrix} \cos(\theta\zeta) \\ \frac{\lambda - \theta^2}{\theta}\sin(\theta\zeta) \end{bmatrix}$
- **虚部**：$\begin{bmatrix} \sin(\theta\zeta) \\ -\frac{\lambda - \theta^2}{\theta}\cos(\theta\zeta) \end{bmatrix}$

共轭根 $m = -\theta i$ 给出相同的实部和虚部（线性相关），因此 $\pm\theta i$ 贡献两个线性无关的实解。

#### 5.4.3 $\pm\omega i$ 根的振型分量

$\pm\omega i$ 根的推导与 Region I 中完全相同（见 5.2.2 节），得到：
- **实部**：$\begin{bmatrix} \cos(\omega\zeta) \\ \frac{\lambda - \omega^2}{\omega}\sin(\omega\zeta) \end{bmatrix}$
- **虚部**：$\begin{bmatrix} \sin(\omega\zeta) \\ -\frac{\lambda - \omega^2}{\omega}\cos(\omega\zeta) \end{bmatrix}$

#### 5.4.4 Region III 的通解

四个纯虚根产生四个线性无关的实解，通解为：

$$
\boxed{\begin{bmatrix} U^*(\zeta) \\ \Phi^*(\zeta) \end{bmatrix} = A_1\begin{bmatrix} \sin\theta\zeta \\ -\frac{\lambda-\theta^2}{\theta}\cos\theta\zeta \end{bmatrix} + A_2\begin{bmatrix} \cos\theta\zeta \\ \frac{\lambda-\theta^2}{\theta}\sin\theta\zeta \end{bmatrix} + A_3\begin{bmatrix} \sin\omega\zeta \\ -\frac{\lambda-\omega^2}{\omega}\cos\omega\zeta \end{bmatrix} + A_4\begin{bmatrix} \cos\omega\zeta \\ \frac{\lambda-\omega^2}{\omega}\sin\omega\zeta \end{bmatrix}} \tag{5.4}
$$

> **注**：Region III 的通解与 Region I 形式上的联系可以通过"解析延拓"理解：将 Region I 中的 $\mu$ 替换为 $i\theta$（因为 $\mu^2 \to -\theta^2$），则 $\sinh(\mu\zeta) \to i\sin(\theta\zeta)$，$\cosh(\mu\zeta) \to \cos(\theta\zeta)$。

---

## 6. 边界条件与特征方程

### 6.1 四种基本边界条件

| 条件名称 | 缩写 | 物理约束 | 方程1 | 方程2 |
|---------|-----|---------|-------|-------|
| Clamped（固定） | C | 位移和转角为零 | $u = 0$ | $\phi = 0$ |
| Pinned（简支） | P | 位移为零，弯矩为零 | $u = 0$ | $\phi' = 0$ |
| Free（自由） | F | 剪力为零，弯矩为零 | $u' - \phi = 0$ | $\phi' = 0$ |
| Roller（滚支） | R | 剪力为零，转角为零 | $u' - \phi = 0$ | $\phi = 0$ |

### 6.2 特征方程推导通用框架

所有边界条件的特征方程推导遵循**统一框架**：

**Step 1**：根据边界类型，在 $\zeta = 0$ 和 $\zeta = 1$ 处写出4个方程

| 边界类型 | $\zeta = 0$ 的方程 | $\zeta = 1$ 的方程 |
|---------|-------------------|-------------------|
| C (固定) | $U^*(0) = 0$, $\Phi^*(0) = 0$ | $U^*(1) = 0$, $\Phi^*(1) = 0$ |
| P (简支) | $U^*(0) = 0$, ${\Phi^*}'(0) = 0$ | $U^*(1) = 0$, ${\Phi^*}'(1) = 0$ |
| F (自由) | ${U^*}'(0) - \Phi^*(0) = 0$, ${\Phi^*}'(0) = 0$ | ${U^*}'(1) - \Phi^*(1) = 0$, ${\Phi^*}'(1) = 0$ |
| R (滚支) | ${U^*}'(0) = 0$, $\Phi^*(0) = 0$ | ${U^*}'(1) = 0$, $\Phi^*(1) = 0$ |

**Step 2**：将通解代入，得到关于 $(A_1, A_2, A_3, A_4)$ 的齐次线性方程组

$$
\mathbf{B} \begin{bmatrix} A_1 \\ A_2 \\ A_3 \\ A_4 \end{bmatrix} = \begin{bmatrix} 0 \\ 0 \\ 0 \\ 0 \end{bmatrix}
$$

**Step 3**：非平凡解存在的条件为 $\det(\mathbf{B}) = 0$

**Step 4**：展开行列式，化简得到特征方程

> **注**：对于不同的区域（Region I/II/III），需要使用不同的通解形式，因此特征方程也会有所不同。

---

### 6.3 Region I ($\lambda < \alpha$) 的特征方程

在此区域，通解为 (5.2) 式，包含双曲函数 $\sinh\mu\zeta$、$\cosh\mu\zeta$ 和三角函数 $\sin\omega\zeta$、$\cos\omega\zeta$。

#### 6.3.1 Pinned-Pinned (PP) 边界

PP边界是最简单的情况，完整推导如下。

**边界条件**：
- 左端 $\zeta = 0$：$U^*(0) = 0$，${\Phi^*}'(0) = 0$
- 右端 $\zeta = 1$：$U^*(1) = 0$，${\Phi^*}'(1) = 0$

**Step 1**：写出通解及其导数

由通解 (5.2)，$U^*$ 和 $\Phi^*$ 及其导数为：
$$
U^*(\zeta) = A_1\sinh\mu\zeta + A_2\cosh\mu\zeta + A_3\sin\omega\zeta + A_4\cos\omega\zeta
$$
$$
\Phi^*(\zeta) = A_1\frac{\lambda+\mu^2}{\mu}\cosh\mu\zeta + A_2\frac{\lambda+\mu^2}{\mu}\sinh\mu\zeta - A_3\frac{\lambda-\omega^2}{\omega}\cos\omega\zeta + A_4\frac{\lambda-\omega^2}{\omega}\sin\omega\zeta
$$
$$
{\Phi^*}'(\zeta) = A_1(\lambda+\mu^2)\sinh\mu\zeta + A_2(\lambda+\mu^2)\cosh\mu\zeta + A_3(\lambda-\omega^2)\sin\omega\zeta + A_4(\lambda-\omega^2)\cos\omega\zeta
$$

**Step 2**：构造边界条件矩阵

将四个边界条件代入通解，得到 $4 \times 4$ 边界条件矩阵 $\mathbf{B}$：

$$
\mathbf{B} = \begin{bmatrix}
0 & 1 & 0 & 1 \\
0 & \lambda+\mu^2 & 0 & \lambda-\omega^2 \\
\sinh\mu & \cosh\mu & \sin\omega & \cos\omega \\
(\lambda+\mu^2)\sinh\mu & (\lambda+\mu^2)\cosh\mu & (\lambda-\omega^2)\sin\omega & (\lambda-\omega^2)\cos\omega
\end{bmatrix}
$$

**Step 3**：计算行列式

按第一列展开，利用余子式计算 $\det(\mathbf{B}) = 0$：

$$
\det(\mathbf{B}) = -\sinh\mu \cdot M_{31} + (\lambda+\mu^2)\sinh\mu \cdot M_{41}
$$

其中 $M_{31}$、$M_{41}$ 为对应的 $3 \times 3$ 余子式。经过完整展开和化简，得到：

$$
\det(\mathbf{B}) = (-\sin\omega)(\mu^2+\omega^2)^2 \sinh\mu = 0
$$

**传统形式（参考文献 Table 7）**：
$$
\boxed{(-\sin\omega)(\mu^2+\omega^2)^2 \sinh\mu = 0} \tag{6.1}
$$

**简化分析**：由于 $(\mu^2+\omega^2)^2 > 0$ 恒成立，且 $\sinh\mu \neq 0$（当 $\mu > 0$，即 $\lambda < \alpha$），故特征方程等价于：
$$
\sin\omega = 0
$$

**特征值**：$\omega_n = n\pi$，$n = 1, 2, 3, \ldots$

**Step 4**：求解振型系数

将特征值代回边界条件，求系数 $(A_1, A_2, A_3, A_4)$。

由 $U^*(0) = 0$：
$$
A_2 + A_4 = 0 \quad \Rightarrow \quad A_4 = -A_2
$$

由 ${\Phi^*}'(0) = 0$：
$$
A_2(\lambda+\mu^2) + A_4(\lambda-\omega^2) = 0
$$
代入 $A_4 = -A_2$：
$$
A_2[(\lambda+\mu^2) - (\lambda-\omega^2)] = A_2(\mu^2+\omega^2) = 0
$$
由于 $\mu^2+\omega^2 > 0$，故 $A_2 = 0$，从而 $A_4 = 0$。

由 $U^*(1) = 0$（此时 $\sin\omega = 0$ 已满足）：
$$
A_1\sinh\mu + A_3\sin\omega = A_1\sinh\mu = 0
$$
由于 $\sinh\mu \neq 0$，故 $A_1 = 0$。

剩余自由度 $A_3$ 用于归一化，取 $A_3 = 1$。

**PP边界的模态振型**：
$$
\boxed{U^*(\zeta) = \sin(\omega_n\zeta), \quad \Phi^*(\zeta) = -\frac{\lambda-\omega_n^2}{\omega_n}\cos(\omega_n\zeta)} \tag{6.2}
$$

> **注**：其他边界条件的振型系数推导方法相同，但系数非零，且传统形式含 $\sinh\mu$、$\cosh\mu$ 会数值溢出。稳定化的振型系数见§7。

#### 6.3.2 Region I 特征方程汇总

| 边界条件 | 传统形式特征方程 | 公式 |
|---------|-----------------|------|
| PP | $(-\sin\omega)(\mu^2+\omega^2)^2 \sinh\mu = 0$ | (6.1) |
| CC | $\sinh\mu\sin\omega\left[\frac{\omega(\lambda+\mu^2)}{\mu(\lambda-\omega^2)} - \frac{\mu(\lambda-\omega^2)}{\omega(\lambda+\mu^2)}\right] + 2\cosh\mu\cos\omega - 2 = 0$ | (6.3) |
| CF | $\left[\frac{\lambda+\mu^2}{\lambda-\omega^2} + \frac{\lambda-\omega^2}{\lambda+\mu^2}\right]\cosh\mu\cos\omega + \left[\frac{\omega}{\mu} - \frac{\mu}{\omega}\right]\sinh\mu\sin\omega - 2 = 0$ | (6.4) |
| CP | $\frac{\mu(\lambda-\omega^2)}{\omega(\lambda+\mu^2)}\sinh\mu\cos\omega + \cosh\mu\sin\omega = 0$ | (6.5) |
| FF | $-\left[\frac{\mu(\lambda+\mu^2)}{\omega(\lambda-\omega^2)} - \frac{\omega(\lambda-\omega^2)}{\mu(\lambda+\mu^2)}\right]\sinh\mu\sin\omega + 2\cosh\mu\cos\omega - 2 = 0$ | (6.6) |
| PF | $\lambda\mu(\lambda+\mu^2)\sinh\mu\cos\omega + \lambda\omega(\lambda-\omega^2)\cosh\mu\sin\omega = 0$ | (6.7) |
| PR | $\cos\omega\cosh\mu = 0$ | (6.8) |
| CR | $\sin\omega\cosh\mu - \frac{\omega(\lambda+\mu^2)}{\mu(\lambda-\omega^2)}\cos\omega\sinh\mu = 0$ | (6.9) |
| RF | $\sin\omega\cosh\mu - \frac{\omega(\lambda-\omega^2)}{\mu(\lambda+\mu^2)}\cos\omega\sinh\mu = 0$ | (6.10) |
| RR | $\sin\omega\sinh\mu = 0$ | (6.11) |

> **Roller边界说明**：Roller（滚支）约束截面转角（$\phi = 0$）但允许横向滑动（$u' = 0$）。工程应用：对称结构的对称面。

---

### 6.4 Region II ($\lambda = \alpha$) 的特征方程

在临界频率 $\lambda = \alpha$ 处，$\Delta = 1$，$\mu = 0$（二重根）。

#### 6.4.1 Region II 通解回顾

由 §5.3，Region II 通解为：
$$
U^*(\zeta) = A_2 + A_3\sin\omega\zeta + A_4\cos\omega\zeta
$$
$$
\Phi^*(\zeta) = A_1 + A_2\alpha\zeta - A_3\frac{\lambda-\omega^2}{\omega}\cos\omega\zeta + A_4\frac{\lambda-\omega^2}{\omega}\sin\omega\zeta
$$

导数：
$$
{U^*}'(\zeta) = A_3\omega\cos\omega\zeta - A_4\omega\sin\omega\zeta
$$
$$
{\Phi^*}'(\zeta) = A_2\alpha + A_3(\lambda-\omega^2)\sin\omega\zeta + A_4(\lambda-\omega^2)\cos\omega\zeta
$$

#### 6.4.2 CC边界推导示例

**边界条件**：$U^*(0)=0$, $\Phi^*(0)=0$, $U^*(1)=0$, $\Phi^*(1)=0$

由 $U^*(0) = A_2 + A_4 = 0$ 得 $A_4 = -A_2$

由 $\Phi^*(0) = A_1 - A_3\frac{\lambda-\omega^2}{\omega} = 0$ 得 $A_1 = A_3\frac{\lambda-\omega^2}{\omega}$

由 $U^*(1) = A_2 + A_3\sin\omega + A_4\cos\omega = 0$，代入 $A_4 = -A_2$：
$$
A_2(1 - \cos\omega) + A_3\sin\omega = 0
$$

由 $\Phi^*(1) = A_1 + A_2\alpha - A_3\frac{\lambda-\omega^2}{\omega}\cos\omega + A_4\frac{\lambda-\omega^2}{\omega}\sin\omega = 0$

代入 $A_1$、$A_4$ 的表达式并整理，联立方程组得到：
$$
\det\begin{bmatrix} 1-\cos\omega & \sin\omega \\ \alpha & \frac{\lambda-\omega^2}{\omega}(1-\cos\omega) \end{bmatrix} = 0
$$

展开：$(1-\cos\omega)\frac{\lambda-\omega^2}{\omega}(1-\cos\omega) - \alpha\sin\omega = 0$

由于 $\lambda = \alpha$，在特征值处 $\frac{\lambda-\omega^2}{\omega} \neq 0$（一般情况），特征方程简化为：
$$
\boxed{\sin\omega(1 - \cos\omega) = 0}
$$

#### 6.4.3 Region II 特征方程汇总

| 边界 | 特征方程 | 边界 | 特征方程 |
|-----|---------|-----|---------|
| PP | $\sin\omega = 0$ | PR | $\cos\omega = 0$ |
| CC | $\sin\omega(1 - \cos\omega) = 0$ | CR | $\sin\omega = 0$ |
| CF | $\cos^2\omega - 1 = 0$ | RF | $\sin\omega = 0$ |
| CP | $\sin\omega = 0$ | RR | $\sin\omega = 0$ |
| FF | $\sin\omega(1 - \cos\omega) = 0$ | | |
| PF | $\cos\omega\sin\omega = 0$ | | |

> **注**：Region II 是孤立点（$\lambda = \alpha$），数值计算中通常跳过此点。

---

### 6.5 Region III ($\lambda > \alpha$) 的特征方程

#### 6.5.1 从 Region I 到 Region III 的解析延拓

**为什么需要变换？**

回顾§4.5中 $\mu^2$ 的定义：
$$
\mu^2 = \frac{1}{2}\lambda(1+\gamma)(\sqrt{\Delta} - 1)
$$

- 当 $\lambda < \alpha$ 时，$\Delta > 1$，故 $\sqrt{\Delta} > 1$，$\mu^2 > 0$，$\mu$ 为**实数**
- 当 $\lambda > \alpha$ 时，$\Delta < 1$，故 $\sqrt{\Delta} < 1$，$\mu^2 < 0$，$\mu$ 为**纯虚数**

设 $\mu = i\theta$，其中 $\theta$ 为实数，则：
$$
\theta^2 = -\mu^2 = \frac{1}{2}\lambda(1+\gamma)(1 - \sqrt{\Delta}) > 0
$$

这正是§4.5中定义的 (4.8)。

**双曲函数到三角函数的变换**：

由复变函数的恒等式：
$$
\sinh(i\theta) = i\sin\theta, \quad \cosh(i\theta) = \cos\theta
$$

因此，在 Region I 的特征方程中：
- $\sinh\mu \to \sinh(i\theta) = i\sin\theta$
- $\cosh\mu \to \cosh(i\theta) = \cos\theta$
- $\mu^2 \to (i\theta)^2 = -\theta^2$，故 $\lambda + \mu^2 \to \lambda - \theta^2$

#### 6.5.2 变换规则汇总

| Region I | Region III | 说明 |
|----------|------------|------|
| $\mu$ | $i\theta$ | 实根变为虚根 |
| $\mu^2$ | $-\theta^2$ | |
| $\sinh\mu$ | $i\sin\theta$ | 双曲正弦 → 虚单位×正弦 |
| $\cosh\mu$ | $\cos\theta$ | 双曲余弦 → 余弦 |
| $\lambda + \mu^2$ | $\lambda - \theta^2$ | 分母中的组合项 |

#### 6.5.3 CC边界推导示例

以 Clamped-Clamped 边界为例，展示从 Region I 到 Region III 的完整变换过程。

**Region I 的 CC 特征方程**（6.3式）：
$$
\sinh\mu\sin\omega\left[\frac{\omega(\lambda+\mu^2)}{\mu(\lambda-\omega^2)} - \frac{\mu(\lambda-\omega^2)}{\omega(\lambda+\mu^2)}\right] + 2\cosh\mu\cos\omega - 2 = 0
$$

**Step 1**：应用变换规则 $\mu \to i\theta$

$$
\sinh(i\theta)\sin\omega\left[\frac{\omega(\lambda+(i\theta)^2)}{i\theta(\lambda-\omega^2)} - \frac{i\theta(\lambda-\omega^2)}{\omega(\lambda+(i\theta)^2)}\right] + 2\cosh(i\theta)\cos\omega - 2 = 0
$$

**Step 2**：代入恒等式 $\sinh(i\theta) = i\sin\theta$，$\cosh(i\theta) = \cos\theta$，$(i\theta)^2 = -\theta^2$

$$
i\sin\theta\sin\omega\left[\frac{\omega(\lambda-\theta^2)}{i\theta(\lambda-\omega^2)} - \frac{i\theta(\lambda-\omega^2)}{\omega(\lambda-\theta^2)}\right] + 2\cos\theta\cos\omega - 2 = 0
$$

**Step 3**：简化方括号内的项

第一项：$\frac{\omega(\lambda-\theta^2)}{i\theta(\lambda-\omega^2)} = -\frac{i\omega(\lambda-\theta^2)}{\theta(\lambda-\omega^2)}$

第二项：$\frac{i\theta(\lambda-\omega^2)}{\omega(\lambda-\theta^2)}$（保持不变）

代入：
$$
i\sin\theta\sin\omega\left[-\frac{i\omega(\lambda-\theta^2)}{\theta(\lambda-\omega^2)} - \frac{i\theta(\lambda-\omega^2)}{\omega(\lambda-\theta^2)}\right] + 2\cos\theta\cos\omega - 2 = 0
$$

**Step 4**：提取公因子 $-i$

$$
i\sin\theta\sin\omega \cdot (-i)\left[\frac{\omega(\lambda-\theta^2)}{\theta(\lambda-\omega^2)} + \frac{\theta(\lambda-\omega^2)}{\omega(\lambda-\theta^2)}\right] + 2\cos\theta\cos\omega - 2 = 0
$$

由于 $i \cdot (-i) = 1$：
$$
\sin\theta\sin\omega\left[\frac{\omega(\lambda-\theta^2)}{\theta(\lambda-\omega^2)} + \frac{\theta(\lambda-\omega^2)}{\omega(\lambda-\theta^2)}\right] + 2\cos\theta\cos\omega - 2 = 0
$$

**Step 5**：整理为标准形式

将方括号内的和写为负号形式（与论文Table 8一致）：
$$
\boxed{\sin\theta\sin\omega\left[-\frac{\omega(\lambda-\theta^2)}{\theta(\lambda-\omega^2)} - \frac{\theta(\lambda-\omega^2)}{\omega(\lambda-\theta^2)}\right] - 2\cos\theta\cos\omega + 2 = 0}
$$

> **注**：上式与 Step 4 结果等价，只是整体乘以 $-1$ 后移项。论文采用此形式是为了与数值稳定化时的符号约定一致。

#### 6.5.4 Region III 特征方程汇总

按照变换规则，得到所有边界条件的 Region III 特征方程：

| 边界 | 特征方程（传统形式） |
|-----|---------|
| PP | $(-\sin\omega)(\omega^2-\theta^2)^2\sin\theta = 0$ |
| CC | $\sin\theta\sin\omega\left[-\frac{\omega(\lambda-\theta^2)}{\theta(\lambda-\omega^2)} - \frac{\theta(\lambda-\omega^2)}{\omega(\lambda-\theta^2)}\right] - 2\cos\theta\cos\omega + 2 = 0$ |
| CF | $\left[\frac{\theta^2-\lambda}{\lambda-\omega^2} + \frac{\lambda-\omega^2}{\theta^2-\lambda}\right]\cos\theta\cos\omega - \frac{\theta^2+\omega^2}{\theta\omega}\sin\theta\sin\omega + 2 = 0$ |
| CP | $(\omega^2-\theta^2)\left[\sin\theta\cos\omega - \frac{\omega(\lambda-\theta^2)}{\theta(\lambda-\omega^2)}\cos\theta\sin\omega\right] = 0$ |
| FF | $\sin\theta\sin\omega\left[-\frac{\theta(\lambda-\theta^2)}{2\omega(\lambda-\omega^2)} - \frac{\omega(\lambda-\omega^2)}{2\theta(\lambda-\theta^2)}\right] - \cos\theta\cos\omega + 1 = 0$ |
| PF | $\frac{\lambda-\omega^2}{\theta}\cos\theta\sin\omega - \frac{\lambda-\theta^2}{\omega}\cos\omega\sin\theta = 0$ |
| PR | $\cos\omega\cos\theta = 0$ |
| CR | $\sin\omega\cos\theta - \frac{\omega(\lambda-\theta^2)}{\theta(\lambda-\omega^2)}\cos\omega\sin\theta = 0$ |
| RF | $\sin\omega\cos\theta - \frac{\omega(\lambda-\omega^2)}{\theta(\lambda-\theta^2)}\cos\omega\sin\theta = 0$ |
| RR | $\sin\omega\sin\theta = 0$ |

> **PP边界的重要说明**（参见参考文献附录B）：
>
> 与 Region I 不同，Region III 的 PP 边界特征方程 $(-\sin\omega)(\omega^2-\theta^2)^2\sin\theta = 0$ 有**两组独立的特征值分支**：
>
> 1. **分支1**：$\sin\omega = 0$（即 $\omega_n = n\pi$）
>    - 振型：$U^*(\zeta) = \sin(\theta\zeta)$，$\Phi^*(\zeta) = -\frac{\lambda-\theta^2}{\theta}\cos(\theta\zeta)$
>
> 2. **分支2**：$\sin\theta = 0$（即 $\theta_n = n\pi$）
>    - 振型：$U^*(\zeta) = \sin(\omega\zeta)$，$\Phi^*(\zeta) = -\frac{\lambda-\omega^2}{\omega}\cos(\omega\zeta)$
>
> 这体现了 Timoshenko 梁在高频区域的**第二频谱**现象。类似地，RR、PR 等边界也具有双分支结构。

---

## 7. 数值稳定化

### 7.1 问题描述

传统特征方程（如 6.3-6.7）在高阶模态时存在严重的**数值不稳定性**。

**原因**：当 $\mu$ 较大时：
- $\sinh(\mu) \approx \cosh(\mu) \approx \frac{1}{2}e^{\mu}$ 呈指数增长
- 不同项之间出现大数相减，导致**灾难性抵消**（catastrophic cancellation）
- 计算结果出现NaN或严重精度损失

**示例**：对于 CC 边界的第25阶模态，$\mu \approx 25$，$\sinh(\mu) \approx 10^{10}$，而最终结果应为0，数值计算无法实现这一精度。

### 7.2 稳定化恒等式

定义以下辅助函数：

**衰减函数 $f(\mu)$**：
$$
\boxed{f(\mu) = \frac{1}{\cosh(\mu)} = \frac{2e^{-\mu}}{1+e^{-2\mu}} \to 0 \quad \text{as } \mu \to \infty} \tag{7.1}
$$

**归一化因子 $g(\mu)$**：
$$
\boxed{g(\mu) = \frac{1}{\sqrt{\cosh^2(\mu)+\sinh^2(\mu)}} = \frac{\sqrt{2}e^{-\mu}}{\sqrt{1+e^{-4\mu}}} \to 0 \quad \text{as } \mu \to \infty} \tag{7.2}
$$

**稳定化角度 $\psi(\mu)$**：
$$
\boxed{\psi(\mu) = \text{atan2}(1-e^{-2\mu}, 1+e^{-2\mu}) \to \frac{\pi}{4} \quad \text{as } \mu \to \infty} \tag{7.3}
$$

### 7.3 关键恒等式推导

**Step 1**：由双曲函数定义
$$
\sinh(\mu) = \frac{e^{\mu} - e^{-\mu}}{2}, \quad \cosh(\mu) = \frac{e^{\mu} + e^{-\mu}}{2}
$$

**Step 2**：计算归一化比值
$$
\frac{\sinh(\mu)}{\sqrt{\cosh^2(\mu)+\sinh^2(\mu)}} = \frac{e^{\mu} - e^{-\mu}}{\sqrt{(e^{\mu}+e^{-\mu})^2 + (e^{\mu}-e^{-\mu})^2}}
$$

分子分母同除以 $e^{\mu}$：
$$
= \frac{1 - e^{-2\mu}}{\sqrt{(1+e^{-2\mu})^2 + (1-e^{-2\mu})^2}}
$$

由于 $(1+a)^2 + (1-a)^2 = 2(1+a^2)$：
$$
= \frac{1 - e^{-2\mu}}{\sqrt{2(1+e^{-4\mu})}} = \frac{\sqrt{2}(1 - e^{-2\mu})}{2\sqrt{1+e^{-4\mu}}}
$$

**Step 3**：引入 $\psi$ 的定义

由 $\psi = \text{atan2}(1-e^{-2\mu}, 1+e^{-2\mu})$，有：
$$
\sin(\psi) = \frac{1-e^{-2\mu}}{\sqrt{(1-e^{-2\mu})^2+(1+e^{-2\mu})^2}} = \frac{1-e^{-2\mu}}{\sqrt{2(1+e^{-4\mu})}}
$$
$$
\cos(\psi) = \frac{1+e^{-2\mu}}{\sqrt{2(1+e^{-4\mu})}}
$$

因此：
$$
\boxed{\frac{\sinh(\mu)}{\sqrt{\cosh^2(\mu)+\sinh^2(\mu)}} = \sin(\psi(\mu))} \tag{7.4}
$$

$$
\boxed{\frac{\cosh(\mu)}{\sqrt{\cosh^2(\mu)+\sinh^2(\mu)}} = \cos(\psi(\mu))} \tag{7.5}
$$

### 7.4 CC边界稳定化推导示例

以 Clamped-Clamped 边界为例，详细展示特征方程和振型系数的稳定化过程。

#### 7.4.1 特征方程的稳定化

**传统形式**（6.3式）：
$$
\sinh\mu\sin\omega\left[\frac{\omega(\lambda+\mu^2)}{\mu(\lambda-\omega^2)} - \frac{\mu(\lambda-\omega^2)}{\omega(\lambda+\mu^2)}\right] + 2\cosh\mu\cos\omega - 2 = 0
$$

**Step 1**：两边同除以 $\sqrt{\cosh^2\mu + \sinh^2\mu}$

由恒等式 (7.4) 和 (7.5)：
$$
\frac{\sinh\mu}{\sqrt{\cosh^2\mu + \sinh^2\mu}} = \sin\psi, \quad \frac{\cosh\mu}{\sqrt{\cosh^2\mu + \sinh^2\mu}} = \cos\psi
$$

代入：
$$
\sin\psi\sin\omega\left[\frac{\omega(\lambda+\mu^2)}{\mu(\lambda-\omega^2)} - \frac{\mu(\lambda-\omega^2)}{\omega(\lambda+\mu^2)}\right] + 2\cos\psi\cos\omega - \frac{2}{\sqrt{\cosh^2\mu + \sinh^2\mu}} = 0
$$

**Step 2**：利用 $g(\mu)$ 定义

由 (7.2) 式：
$$
g(\mu) = \frac{1}{\sqrt{\cosh^2\mu + \sinh^2\mu}} = \frac{\sqrt{2}e^{-\mu}}{\sqrt{1+e^{-4\mu}}}
$$

**Step 3**：得到稳定形式
$$
\boxed{\sin\psi\sin\omega\left[\frac{\omega(\lambda+\mu^2)}{\mu(\lambda-\omega^2)} - \frac{\mu(\lambda-\omega^2)}{\omega(\lambda+\mu^2)}\right] + 2\cos\psi\cos\omega - 2g(\mu) = 0} \tag{7.6}
$$

**数值稳定性分析**：
- $\sin\psi \to \frac{\sqrt{2}}{2}$，$\cos\psi \to \frac{\sqrt{2}}{2}$ 当 $\mu \to \infty$（有界）
- $g(\mu) \to 0$ 当 $\mu \to \infty$（指数衰减）
- 整个表达式有界，不存在大数相减

#### 7.4.2 振型系数的稳定化

**传统形式的振型**（取 $A_1 = 1$）：
$$
U^*(\zeta) = \sinh(\mu\zeta) + A_2\cosh(\mu\zeta) + A_3\sin(\omega\zeta) + A_4\cos(\omega\zeta)
$$

其中 $A_2 = -\left(1 + P(\mu)e^{-\mu}\right)$，$A_4 = R \cdot A_2$。

**问题**：当 $\mu$ 较大时，$\sinh(\mu\zeta)$ 和 $\cosh(\mu\zeta)$ 指数增长，但两者相消后结果有界。

**稳定化策略**：利用恒等式重写双曲函数组合

$$
\sinh(\mu\zeta) + A_2\cosh(\mu\zeta) = -e^{-\mu\zeta} - \frac{P(\mu)}{2}\left(e^{-\mu(1-\zeta)} + e^{-\mu(1+\zeta)}\right)
$$

其中 $P(\mu)$ 的稳定表达式为：
$$
P(\mu) = 2\frac{-e^{-\mu} + \sin\omega\frac{\omega(\lambda+\mu^2)}{\mu(\lambda-\omega^2)} + \cos\omega}{(1+e^{-2\mu}) - 2e^{-\mu}\cos\omega}
$$

**稳定化振型**：
$$
\boxed{U^*(\zeta) = -e^{-\mu\zeta} - \frac{P}{2}\left(e^{-\mu(1-\zeta)} + e^{-\mu(1+\zeta)}\right) + A_3\sin(\omega\zeta) - R(1+Pe^{-\mu})\cos(\omega\zeta)} \tag{7.7}
$$

$$
\boxed{\Phi^*(\zeta) = \frac{\lambda+\mu^2}{\mu}\left[e^{-\mu\zeta} - \frac{P}{2}\left(e^{-\mu(1-\zeta)} - e^{-\mu(1+\zeta)}\right)\right] - A_3\frac{\lambda-\omega^2}{\omega}\cos(\omega\zeta) - R(1+Pe^{-\mu})\frac{\lambda-\omega^2}{\omega}\sin(\omega\zeta)} \tag{7.8}
$$

### 7.5 稳定化公式汇总

#### 7.5.1 稳定形式特征方程汇总（$\lambda < \alpha$）

| 边界 | 稳定形式特征方程 |
|-----|----------------|
| PP | $\sin\omega = 0$（本身已稳定） |
| CC | $\sin\psi\sin\omega\left[\frac{\omega(\lambda+\mu^2)}{\mu(\lambda-\omega^2)} - \frac{\mu(\lambda-\omega^2)}{\omega(\lambda+\mu^2)}\right] + 2\cos\psi\cos\omega - 2g = 0$ |
| CF | $\left[\frac{\lambda+\mu^2}{\lambda-\omega^2} + \frac{\lambda-\omega^2}{\lambda+\mu^2}\right]\cos\psi\cos\omega + \left[\frac{\omega}{\mu} - \frac{\mu}{\omega}\right]\sin\psi\sin\omega - 2g = 0$ |
| CP | $\frac{\mu(\lambda-\omega^2)}{\omega(\lambda+\mu^2)}\sin\psi\cos\omega + \cos\psi\sin\omega = 0$ |
| FF | $-\left[\frac{\mu(\lambda+\mu^2)}{\omega(\lambda-\omega^2)} - \frac{\omega(\lambda-\omega^2)}{\mu(\lambda+\mu^2)}\right]\sin\psi\sin\omega + 2\cos\psi\cos\omega - 2g = 0$ |
| PF | $\sqrt{\frac{\mu}{\omega}}\sin\psi\cos\omega + \sqrt{\frac{\omega}{\mu}}\frac{\lambda-\omega^2}{\lambda+\mu^2}\cos\psi\sin\omega = 0$ |

| 边界 | 稳定形式特征方程 |
|-----|----------------|
| PR | $\cos\omega\cos\psi = 0$ |
| CR | $\sin\omega\cos\psi - \frac{\omega(\lambda+\mu^2)}{\mu(\lambda-\omega^2)}\cos\omega\sin\psi = 0$ |
| RF | $\sin\omega\cos\psi - \frac{\omega(\lambda-\omega^2)}{\mu(\lambda+\mu^2)}\cos\omega\sin\psi = 0$ |
| RR | $\sin\omega\sin\psi = 0$ |

#### 7.5.2 稳定化振型系数汇总（$\lambda < \alpha$）

振型通式为（取 $A_1 = 1$ 归一化）：
$$
U^*(\zeta) = -e^{-\mu\zeta} - \frac{P}{2}\left(e^{-\mu(1-\zeta)} + e^{-\mu(1+\zeta)}\right) + A_3\sin(\omega\zeta) + R(1+Pe^{-\mu})\cos(\omega\zeta)
$$

| 边界 | $P(\mu, \omega, \lambda)$ | $A_3$ | $R$ |
|-----|--------------------------|-------|-----|
| CC | $2\frac{-e^{-\mu}+\sin\omega\frac{\omega(\lambda+\mu^2)}{\mu(\lambda-\omega^2)}+\cos\omega}{(1+e^{-2\mu})-2e^{-\mu}\cos\omega}$ | $\frac{\omega(\mu^2+\lambda)}{\mu(\lambda-\omega^2)}$ | $-1$ |
| CF | $2\frac{-\omega e^{-\mu}(\lambda-\omega^2)+\omega(\lambda+\mu^2)\cos\omega-\mu(\lambda-\omega^2)\sin\omega}{(\lambda-\omega^2)(\omega(e^{-2\mu}-1)+2\mu e^{-\mu}\sin\omega)}$ | $\frac{\omega(\mu^2+\lambda)}{\mu(\lambda-\omega^2)}$ | $-1$ |
| CP | $2\frac{-\mu e^{-\mu}(\lambda-\omega^2)+\omega(\lambda+\mu^2)\sin\omega+\mu(\lambda-\omega^2)\cos\omega}{\mu(e^{-2\mu}+1)(\lambda-\omega^2)-2\mu e^{-\mu}(\lambda-\omega^2)\cos\omega}$ | $\frac{\omega(\mu^2+\lambda)}{\mu(\lambda-\omega^2)}$ | $-1$ |
| FF | $-2\frac{\omega e^{-\mu}(\lambda-\omega^2)+\mu(\lambda+\mu^2)\sin\omega-\omega(\lambda-\omega^2)\cos\omega}{2\mu e^{-\mu}(\lambda+\mu^2)\sin\omega+\omega(e^{-2\mu}-1)(\lambda-\omega^2)}$ | $\frac{\omega}{\mu}$ | $-\frac{\mu^2+\lambda}{\lambda-\omega^2}$ |
| CR | $\frac{(2\omega\mu^2+2\lambda\omega)\cos\omega+(2\mu\omega^2-2\lambda\mu)\sin\omega-2\lambda\omega e^{-\mu}-2\mu^2\omega e^{-\mu}}{-\lambda\omega-\mu^2\omega+(2\lambda\mu e^{-\mu}-2\mu\omega^2 e^{-\mu})\sin\omega+\lambda\omega e^{-2\mu}+\mu^2\omega e^{-2\mu}}$ | $\frac{\omega(\mu^2+\lambda)}{\mu(\lambda-\omega^2)}$ | $-1$ |

**PP、PF、PR、RR边界**（无需稳定化）：

| 边界 | 振型表达式 |
|-----|---------|
| PP | $U^*(\zeta) = \sin(\omega\zeta)$，$\Phi^*(\zeta) = -\frac{\lambda-\omega^2}{\omega}\cos(\omega\zeta)$ |
| PF | $U^*(\zeta) = \frac{\mu\cos\omega}{\omega}\frac{e^{-\mu(1-\zeta)}-e^{-\mu(1+\zeta)}}{1+e^{-2\mu}} + \sin(\omega\zeta)$，$\Phi^*(\zeta) = \frac{(\lambda+\mu^2)\cos\omega}{\omega}\frac{e^{-\mu(1-\zeta)}+e^{-\mu(1+\zeta)}}{1+e^{-2\mu}} - \frac{\lambda-\omega^2}{\omega}\cos(\omega\zeta)$ |
| PR | $U^*(\zeta) = \sin(\omega\zeta)$，$\Phi^*(\zeta) = -\frac{\lambda-\omega^2}{\omega}\cos(\omega\zeta)$ |
| RR | $U^*(\zeta) = \cos(\omega\zeta)$，$\Phi^*(\zeta) = \frac{\lambda-\omega^2}{\omega}\sin(\omega\zeta)$ |

**RF边界**（稳定化形式）：

$$
U^*(\zeta) = \frac{\cos\omega(\omega^2-\lambda)}{(\mu^2+\lambda)}\frac{e^{-\mu(1-\zeta)}+e^{-\mu(1+\zeta)}}{1+e^{-2\mu}} + \cos(\omega\zeta)
$$
$$
\Phi^*(\zeta) = \frac{\cos\omega(\omega^2-\lambda)}{\mu}\frac{e^{-\mu(1-\zeta)}-e^{-\mu(1+\zeta)}}{1+e^{-2\mu}} + \frac{\lambda-\omega^2}{\omega}\sin(\omega\zeta)
$$

### 7.6 稳定化效果对比

| 形式 | 特点 | 适用范围 |
|-----|------|---------|
| **传统形式** | 含 $\sinh$, $\cosh$，指数增长 | 低阶模态 ($\mu < 10$) |
| **稳定形式** | 含 $\sin\psi$, $\cos\psi$，有界振荡 | 任意阶模态 |

稳定形式的特征函数值域有界（大约在 $[-2, 2]$ 范围内），使得数值求根算法能够可靠工作。

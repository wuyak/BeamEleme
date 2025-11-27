# 压缩与扭转模态理论推导

> **参考文献**:
> - Blevins, R. D. (1979). *Formulas for Natural Frequency and Mode Shape*. Van Nostrand Reinhold. ISBN: 0-442-20710-7
>   - 轴向振动：Table 8-16
>   - 扭转振动：Table 8-19

---

## 目录

1. [引言](#1-引言)
2. [波动方程推导](#2-波动方程推导)
3. [分离变量法求解](#3-分离变量法求解)
4. [边界条件与特征方程](#4-边界条件与特征方程)
5. [模态振型与固有频率](#5-模态振型与固有频率)
6. [有限元离散化](#6-有限元离散化)
7. [刚度矩阵与质量矩阵](#7-刚度矩阵与质量矩阵)
8. [统一性总结](#8-统一性总结)

---

## 1. 引言

压缩（轴向）振动和扭转振动是梁振动中最简单的两类模态。与弯曲-剪切耦合模态不同，它们的控制方程都是**标准一维波动方程**，数学结构完全相同，只是物理参数不同。

| 振动类型 | 刚度参数 | 惯性参数 | 自由度 |
|---------|---------|---------|-------|
| 轴向（压缩） | $EA$ | $\rho A$ | 轴向位移 $u$ |
| 扭转 | $GJ$ | $\rho I_p$ | 扭转角 $\phi$ |

**符号说明**：
- $E$：弹性模量 (Pa)
- $G$：剪切模量 (Pa)
- $A$：截面面积 (m²)
- $J$：扭转常数（圆截面等于极惯性矩 $I_p$，非圆截面需另行计算）(m⁴)
- $I_p$：极惯性矩 (m⁴)
- $\rho$：材料密度 (kg/m³)
- $L$：梁长 (m)

---

## 2. 波动方程推导

### 2.1 轴向振动

考虑均匀杆的轴向自由振动，设 $u(x,t)$ 为轴向位移。

**微元受力分析**：取长度为 $dx$ 的微元，设截面轴力为 $N(x,t)$。

根据胡克定律，轴向应力与应变的关系为：
$$
\sigma = E\varepsilon = E\frac{\partial u}{\partial x}
$$

轴力为：
$$
N = \sigma A = EA\frac{\partial u}{\partial x} \tag{2.1}
$$

微元 $[x, x+dx]$ 左端的轴力为 $N(x)$，右端的轴力为 $N(x+dx)$，沿梁长的轴力净增量为：
$$
\Delta N = N(x+dx) - N(x)
$$

微元的质量为：
$$
dm = \rho \cdot A \cdot dx
$$

微元的加速度为 $a = \frac{\partial^2 u}{\partial t^2}$。

由牛顿第二定律 $\sum F = dm \cdot a$：
$$
N(x+dx) - N(x) = \rho A \, dx \cdot \frac{\partial^2 u}{\partial t^2}
$$

$$
\frac{\partial N}{\partial x} dx = \rho A \, dx \cdot \frac{\partial^2 u}{\partial t^2}
$$

将 (2.1) 代入：
$$
\frac{\partial}{\partial x}\left(EA\frac{\partial u}{\partial x}\right) = \rho A \frac{\partial^2 u}{\partial t^2}
$$

对于均匀杆（$E$、$A$ 为常数）：

$$
\boxed{c^2 \frac{\partial^2 u}{\partial x^2} = \frac{\partial^2 u}{\partial t^2}, \quad c = \sqrt{\frac{E}{\rho}}} \tag{2.2}
$$

其中 $c$ 为轴向波速（纵波波速）。

### 2.2 扭转振动

考虑均匀杆的扭转自由振动，设 $\phi(x,t)$ 为扭转角。

**微元力矩分析**：取长度为 $dx$ 的微元，设截面扭矩为 $T(x,t)$。

根据扭转本构关系，扭矩与扭转率的关系为：
$$
T = GJ\frac{\partial \phi}{\partial x} \tag{2.3}
$$

其中 $J$ 为截面抗扭刚度（对于圆截面 $J = I_p$）。

微元 $[x, x+dx]$ 左端的扭矩为 $T(x)$，右端的扭矩为 $T(x+dx)$，沿梁长的扭矩净增量为：
$$
\Delta T = T(x+dx) - T(x)
$$

微元绕轴线的转动惯量为：
$$
dI = \rho \cdot I_p \cdot dx
$$

其中 $I_p$ 为截面极惯性矩。微元的角加速度为 $\alpha = \frac{\partial^2 \phi}{\partial t^2}$。

由角动量定理 $\sum M = dI \cdot \alpha$：
$$
T(x+dx) - T(x) = \rho I_p \, dx \cdot \frac{\partial^2 \phi}{\partial t^2}
$$

$$
\frac{\partial T}{\partial x} dx = \rho I_p \, dx \cdot \frac{\partial^2 \phi}{\partial t^2}
$$

将 (2.3) 代入：
$$
\frac{\partial}{\partial x}\left(GJ\frac{\partial \phi}{\partial x}\right) = \rho I_p \frac{\partial^2 \phi}{\partial t^2}
$$

对于均匀杆：

$$
\boxed{c^2 \frac{\partial^2 \phi}{\partial x^2} = \frac{\partial^2 \phi}{\partial t^2}, \quad c = \sqrt{\frac{GJ}{\rho I_p}}} \tag{2.4}
$$

其中 $c$ 为扭转波速。

### 2.3 统一形式

两者具有相同的数学结构——**一维波动方程**：

$$
c^2 \frac{\partial^2 u}{\partial x^2} = \frac{\partial^2 u}{\partial t^2} \tag{2.5}
$$

下面的推导对两种振动都适用，只需将 $u$ 理解为广义位移（轴向位移或扭转角），$c$ 理解为相应的波速。

---

## 3. 分离变量法求解

### 3.1 分离变量假设

假设解的形式为：

$$
u(x,t) = U(x) T(t) \tag{3.1}
$$

代入波动方程 (2.5)：

$$
c^2 U''(x) T(t) = U(x) \ddot{T}(t)
$$

分离变量：

$$
\frac{U''(x)}{U(x)} = \frac{1}{c^2}\frac{\ddot{T}(t)}{T(t)} = -k^2 \tag{3.2}
$$

等式左边只是 $x$ 的函数，右边只是 $t$ 的函数，两者相等意味着它们都等于同一常数。取为 $-k^2$（负值保证振动解）。

### 3.2 空间方程

$$
U''(x) + k^2 U(x) = 0 \tag{3.3}
$$

通解为：

$$
\boxed{U(x) = A\cos(kx) + B\sin(kx)} \tag{3.4}
$$

其中 $A$、$B$ 由边界条件确定，$k$ 为待定的特征值。

### 3.3 时间方程

$$
\ddot{T}(t) + \omega^2 T(t) = 0, \quad \omega = kc \tag{3.5}
$$

通解为：

$$
T(t) = C\cos(\omega t) + D\sin(\omega t) \tag{3.6}
$$

其中 $\omega = kc$ 为圆频率，$C$、$D$ 由初始条件确定。

---

## 4. 边界条件与特征方程

### 4.1 边界条件类型

对于一维振动问题，每端有两种基本边界条件：

| 类型 | 物理意义 | 数学表达 |
|-----|---------|---------|
| **固定 (C)** | 位移为零 | $U = 0$ |
| **自由 (F)** | 应力/力矩为零 | $U' = 0$ |

### 4.2 CC（两端固定）

**边界条件**：$U(0) = 0$，$U(L) = 0$

**Step 1**：由 $U(0) = 0$
$$
A\cos(0) + B\sin(0) = 0 \quad \Rightarrow \quad A = 0
$$

**Step 2**：由 $U(L) = 0$
$$
B\sin(kL) = 0
$$

要使 $B \neq 0$（非平凡解）：

$$
\boxed{\sin(kL) = 0 \quad \Rightarrow \quad k_i L = i\pi, \quad i = 1,2,3,\ldots} \tag{4.1}
$$

### 4.3 CF（左端固定，右端自由）

**边界条件**：$U(0) = 0$，$U'(L) = 0$

**Step 1**：由 $U(0) = 0$
$$
A = 0
$$

**Step 2**：计算 $U'(x)$
$$
U'(x) = Bk\cos(kx)
$$

**Step 3**：由 $U'(L) = 0$
$$
Bk\cos(kL) = 0
$$

要使 $B \neq 0$：

$$
\boxed{\cos(kL) = 0 \quad \Rightarrow \quad k_i L = \frac{(2i-1)\pi}{2}, \quad i = 1,2,3,\ldots} \tag{4.2}
$$

### 4.4 FC（左端自由，右端固定）

**边界条件**：$U'(0) = 0$，$U(L) = 0$

**Step 1**：计算 $U'(x)$
$$
U'(x) = -Ak\sin(kx) + Bk\cos(kx)
$$

**Step 2**：由 $U'(0) = 0$
$$
Bk = 0 \quad \Rightarrow \quad B = 0
$$

**Step 3**：由 $U(L) = 0$
$$
A\cos(kL) = 0
$$

要使 $A \neq 0$：

$$
\boxed{\cos(kL) = 0 \quad \Rightarrow \quad k_i L = \frac{(2i-1)\pi}{2}, \quad i = 1,2,3,\ldots} \tag{4.3}
$$

### 4.5 FF（两端自由）

**边界条件**：$U'(0) = 0$，$U'(L) = 0$

**Step 1**：由 $U'(0) = 0$
$$
B = 0
$$

**Step 2**：由 $U'(L) = 0$
$$
-Ak\sin(kL) = 0
$$

要使 $A \neq 0$：

$$
\boxed{\sin(kL) = 0 \quad \Rightarrow \quad k_i L = i\pi, \quad i = 1,2,3,\ldots} \tag{4.4}
$$

**注意**：FF 边界条件还存在 $k = 0$ 的刚体模态，对应整体平动（轴向）或整体转动（扭转）。

### 4.6 特征值汇总

| 边界条件 | 特征方程 | 特征值 $\lambda_i = k_i L$ |
|---------|---------|---------------------------|
| CC | $\sin(kL) = 0$ | $i\pi$ |
| CF | $\cos(kL) = 0$ | $(2i-1)\pi/2$ |
| FC | $\cos(kL) = 0$ | $(2i-1)\pi/2$ |
| FF | $\sin(kL) = 0$ | $i\pi$（另有 $\lambda_0 = 0$ 刚体模态） |

---

## 5. 模态振型与固有频率

### 5.1 模态振型函数

将特征值代回通解 (3.4)，得到各边界条件下的模态振型：

| 边界条件 | 模态函数 $U_i(x)$ | 特点 |
|---------|------------------|------|
| CC | $\sin\left(\dfrac{i\pi x}{L}\right)$ | 两端为零 |
| CF | $\sin\left(\dfrac{(2i-1)\pi x}{2L}\right)$ | 左端为零，右端极值 |
| FC | $\cos\left(\dfrac{(2i-1)\pi x}{2L}\right)$ | 左端极值，右端为零 |
| FF | $\cos\left(\dfrac{i\pi x}{L}\right)$ | 两端极值 |

FF 边界条件的刚体模态（$i=0$）：$U_0(x) = 1$（常数）

### 5.2 固有频率公式

由 $\omega_i = k_i c = \frac{\lambda_i}{L} c$，固有频率为：

$$
f_i = \frac{\omega_i}{2\pi} = \frac{\lambda_i}{2\pi L} c \tag{5.1}
$$

**轴向振动**：
$$
\boxed{f_i = \frac{\lambda_i}{2\pi L} \sqrt{\frac{E}{\rho}}} \tag{5.2}
$$

**扭转振动**：
$$
\boxed{f_i = \frac{\lambda_i}{2\pi L} \sqrt{\frac{GJ}{\rho I_p}}} \tag{5.3}
$$

---

## 6. 有限元离散化

上述解析解仅适用于均匀杆的标准边界条件。对于复杂几何或边界条件，需要采用有限元方法。

有限元法的核心思想：
1. 将连续体离散为若干单元
2. 在每个单元内用形函数插值位移场
3. 通过变分原理导出单元刚度矩阵和质量矩阵
4. 组装全局矩阵，求解广义特征值问题

对于轴向和扭转振动，采用最简单的**两节点线性单元**。

---

## 7. 刚度矩阵与质量矩阵

### 7.1 形函数

采用线性形函数，无量纲坐标 $\xi = x/L \in [0,1]$：
$$
N_1(\xi) = 1 - \xi, \quad N_2(\xi) = \xi \tag{7.1}
$$

位移插值：
$$
u(\xi) = N_1 u_1 + N_2 u_2 = [N]\{d\} \tag{7.2}
$$

其中 $[N] = [1-\xi, \xi]$，$\{d\} = \{u_1, u_2\}^T$。

### 7.2 应变-位移矩阵

应变为：
$$
\varepsilon = \frac{\partial u}{\partial x} = \frac{1}{L}\frac{du}{d\xi} = \frac{1}{L}[-1, 1]\{d\} = [B]\{d\} \tag{7.3}
$$

其中 $[B] = \frac{1}{L}[-1, 1]$。

### 7.3 轴向振动

**应变能**：
$$
U = \frac{1}{2}\int_0^L EA \varepsilon^2 \, dx = \frac{1}{2}\{d\}^T \left[\int_0^L [B]^T EA [B] \, dx\right] \{d\} = \frac{1}{2}\{d\}^T [K_e] \{d\}
$$

**刚度矩阵**：
$$
[K_e] = \int_0^L [B]^T EA [B] \, dx = EA \cdot L \int_0^1 \frac{1}{L^2}\begin{bmatrix} 1 \\ -1 \end{bmatrix}[-1, 1] d\xi = \frac{EA}{L}\begin{bmatrix} 1 & -1 \\ -1 & 1 \end{bmatrix}
$$

$$
\boxed{[K_e]_{\text{axial}} = \frac{EA}{L} \begin{bmatrix} 1 & -1 \\ -1 & 1 \end{bmatrix}} \tag{7.4}
$$

**动能**：
$$
T = \frac{1}{2}\int_0^L \rho A \dot{u}^2 \, dx = \frac{1}{2}\{\dot{d}\}^T \left[\int_0^L [N]^T \rho A [N] \, dx\right] \{\dot{d}\} = \frac{1}{2}\{\dot{d}\}^T [M_e] \{\dot{d}\}
$$

**质量矩阵**：
$$
[M_e] = \int_0^L [N]^T \rho A [N] \, dx = \rho A L \int_0^1 \begin{bmatrix} (1-\xi)^2 & (1-\xi)\xi \\ (1-\xi)\xi & \xi^2 \end{bmatrix} d\xi
$$

计算积分：
$$
\int_0^1 (1-\xi)^2 d\xi = \frac{1}{3}, \quad \int_0^1 \xi^2 d\xi = \frac{1}{3}, \quad \int_0^1 (1-\xi)\xi d\xi = \frac{1}{6}
$$

$$
\boxed{[M_e]_{\text{axial}} = \frac{\rho A L}{6} \begin{bmatrix} 2 & 1 \\ 1 & 2 \end{bmatrix}} \tag{7.5}
$$

### 7.4 扭转振动

将轴向振动中的物理量替换：
- $u \to \phi$（广义位移）
- $EA \to GJ$（刚度系数）
- $\rho A \to \rho I_p$（惯性系数）

**刚度矩阵**：
$$
\boxed{[K_e]_{\text{torsion}} = \frac{GJ}{L} \begin{bmatrix} 1 & -1 \\ -1 & 1 \end{bmatrix}} \tag{7.6}
$$

**质量矩阵**：
$$
\boxed{[M_e]_{\text{torsion}} = \frac{\rho I_p L}{6} \begin{bmatrix} 2 & 1 \\ 1 & 2 \end{bmatrix}} \tag{7.7}
$$

---

## 8. 统一性总结

轴向和扭转振动具有完全相同的数学结构：

| | 轴向 | 扭转 |
|---|-----|------|
| 广义位移 | $u$ | $\phi$ |
| 刚度系数 | $EA$ | $GJ$ |
| 惯性系数 | $\rho A$ | $\rho I_p$ |
| 波速 | $\sqrt{E/\rho}$ | $\sqrt{GJ/(\rho I_p)}$ |
| 频率公式 | $\dfrac{\lambda_i}{2\pi L}\sqrt{\dfrac{E}{\rho}}$ | $\dfrac{\lambda_i}{2\pi L}\sqrt{\dfrac{GJ}{\rho I_p}}$ |
| 刚度矩阵 | $\dfrac{EA}{L}\begin{bmatrix} 1 & -1 \\ -1 & 1 \end{bmatrix}$ | $\dfrac{GJ}{L}\begin{bmatrix} 1 & -1 \\ -1 & 1 \end{bmatrix}$ |
| 质量矩阵 | $\dfrac{\rho A L}{6}\begin{bmatrix} 2 & 1 \\ 1 & 2 \end{bmatrix}$ | $\dfrac{\rho I_p L}{6}\begin{bmatrix} 2 & 1 \\ 1 & 2 \end{bmatrix}$ |

这种统一性源于它们都满足一维波动方程 $c^2 u_{xx} = u_{tt}$，只是物理解释不同。

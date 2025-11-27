# Timoshenko梁理论推导

> **核心参考文献**:
> - Khasawneh, F. A., & Segalman, D. J. (2019). Exact and numerically stable expressions for Euler-Bernoulli and Timoshenko beam modes. *International Journal of Mechanical Sciences*, 166, 105234. DOI: https://doi.org/10.1016/j.ijmecsci.2019.105234
> - Friedman, Z., & Kosmatka, J. B. (1993). An improved two-node Timoshenko beam finite element. *Computers & Structures*, 47(3), 473-481. DOI: https://doi.org/10.1016/0045-7949(93)90243-7

---

## 目录

1. [引言与物理背景](#1-引言与物理背景)
2. [Euler-Bernoulli梁理论](#2-euler-bernoulli梁理论)
3. [Timoshenko梁基础](#3-timoshenko梁基础)
4. [Timoshenko梁控制方程](#4-timoshenko梁控制方程)
5. [Timoshenko梁边界条件与模态分析](#5-timoshenko梁边界条件与模态分析)
6. [两种理论的详细对比](#6-两种理论的详细对比)
7. [Hamilton原理推导（能量方法）](#7-hamilton原理推导能量方法)
8. [参考文献](#8-参考文献)

---

## 1. 引言与物理背景

### 1.1 梁振动问题的工程意义

梁是工程结构中最基本的承载构件之一。当梁受到动态荷载或初始扰动时，会产生振动。理解梁的振动特性——包括固有频率和模态振型——对于以下工程问题至关重要：

1. **共振避免**：确保工作频率远离固有频率
2. **结构设计**：优化梁的几何参数以满足动态性能要求
3. **故障诊断**：通过振动响应识别结构损伤
4. **噪声控制**：预测和控制结构振动引起的噪声

### 1.2 两种经典梁理论

历史上，梁的横向振动理论主要有两种：

| 理论 | 提出者 | 年代 | 关键特征 |
|------|--------|------|----------|
| **Euler-Bernoulli理论** | Euler, Bernoulli | 1750 | 忽略剪切变形和转动惯量 |
| **Timoshenko理论** | S.P. Timoshenko | 1921 | 考虑剪切变形和转动惯量 |

**为什么需要Timoshenko理论？**

Euler-Bernoulli理论假设截面变形后仍垂直于中性轴，这意味着忽略了剪切变形。对于细长梁，这个假设是合理的。但对于：
- 短粗梁（细长比$L/h < 10$）
- 高阶振动模态
- 复合材料梁（剪切模量相对较低）
- 夹层结构

剪切变形的影响不可忽略，此时必须采用Timoshenko理论。

### 1.3 坐标系与符号约定

本文采用以下坐标系和符号约定：

**坐标系**：
- $x$轴：沿梁的轴向方向，原点在梁的左端
- $z$轴：垂直于中性轴，向下为正
- 梁长度为$L$，截面高度为$h$

**基本符号**：
| 符号 | 含义 | 单位 |
|------|------|------|
| $w(x,t)$ | 横向位移（挠度） | m |
| $\psi(x,t)$ | 截面转角 | rad |
| $E$ | 弹性模量 | Pa |
| $G$ | 剪切模量 | Pa |
| $\rho$ | 材料密度 | kg/m³ |
| $A$ | 截面面积 | m² |
| $I$ | 截面惯性矩 | m⁴ |
| $\kappa$ | 剪切修正系数 | 无量纲 |
| $L$ | 梁长度 | m |

**内力符号**：
| 符号 | 含义 | 单位 |
|------|------|------|
| $M(x,t)$ | 弯矩 | N·m |
| $V(x,t)$ | 剪力 | N |

**内力正方向约定**：
> 本文采用如下符号约定：
> - **剪力$V$**：对微元左截面，向上为正；对右截面，向下为正
> - **弯矩$M$**：使梁下部纤维受拉为正（即"正弯矩使梁向下凸"）
> - **转角$\psi$**：逆时针转动为正
> - **位移$w$**：向下为正

---

## 2. Euler-Bernoulli梁理论

在推导Timoshenko理论之前，我们首先详细推导Euler-Bernoulli理论，这有助于理解两种理论的本质区别。

### 2.1 基本假设

Euler-Bernoulli梁理论基于以下假设：

**假设1：平截面假设（Plane sections remain plane）**

> 变形前垂直于梁中性轴的平面截面，在变形后仍然保持为平面，且仍然垂直于变形后的中性轴。

这个假设有两层含义：
1. 截面保持平面（不发生翘曲）
2. 截面与中性轴保持垂直（无剪切变形）

**假设2：小变形假设**

> 梁的位移和转角远小于梁的几何尺寸，因此可以采用线性应变-位移关系，且可以在未变形构形上建立平衡方程。

数学表达：
$$
\left| \frac{\partial w}{\partial x} \right| \ll 1, \quad |w| \ll L
$$

**假设3：材料线弹性假设**

> 材料服从广义胡克定律，应力与应变成正比。

$$
\sigma = E \varepsilon
$$

**假设4：忽略轴向惯性力**

> 梁的轴向位移引起的惯性力可以忽略。

**假设5：忽略转动惯量**

> 截面绕其形心转动的动能可以忽略。

### 2.2 运动学关系（位移场）

根据平截面假设，当梁发生弯曲变形时，截面绕中性轴转动。设中性轴上某点的横向位移为$w(x,t)$。

由于截面保持垂直于中性轴，截面的转角$\theta$等于中性轴的斜率：

$$
\theta(x,t) = \frac{\partial w(x,t)}{\partial x} \tag{2.1}
$$

考虑距中性轴距离为$z$的一点，由于截面转动，该点产生轴向位移：

$$
u(x,z,t) = -z \cdot \theta(x,t) = -z \frac{\partial w(x,t)}{\partial x} \tag{2.2}
$$

**物理解释**：负号表示当$\theta > 0$（截面逆时针转动）时，$z > 0$（中性轴下方）的点向$-x$方向移动。

完整的位移场为：
$$
\begin{aligned}
u_x(x,z,t) &= -z \frac{\partial w}{\partial x} \quad &\text{（轴向位移）} \\
u_z(x,z,t) &= w(x,t) \quad &\text{（横向位移）}
\end{aligned} \tag{2.3}
$$

### 2.3 应变-位移关系

由小变形假设，轴向应变为：

$$
\varepsilon_{xx} = \frac{\partial u_x}{\partial x} = \frac{\partial}{\partial x}\left( -z \frac{\partial w}{\partial x} \right) = -z \frac{\partial^2 w}{\partial x^2} \tag{2.4}
$$

定义曲率$\kappa$为：
$$
\kappa = \frac{\partial^2 w}{\partial x^2} \tag{2.5}
$$

则：
$$
\varepsilon_{xx} = -z \kappa \tag{2.6}
$$

**物理解释**：
- 当$\kappa > 0$（向下凸），$z > 0$处（下部）$\varepsilon_{xx} < 0$（受压），$z < 0$处（上部）$\varepsilon_{xx} > 0$（受拉）
- 这与我们对梁弯曲的直观理解一致

剪切应变（在EB理论中）：
$$
\gamma_{xz} = \frac{\partial u_x}{\partial z} + \frac{\partial u_z}{\partial x} = -\frac{\partial w}{\partial x} + \frac{\partial w}{\partial x} = 0 \tag{2.7}
$$

**重要结论**：在Euler-Bernoulli理论中，剪切应变恒为零。这是"截面保持垂直于中性轴"假设的直接数学结果。

### 2.4 本构关系（应力-应变关系）

根据胡克定律，轴向应力为：

$$
\sigma_{xx} = E \varepsilon_{xx} = -Ez \frac{\partial^2 w}{\partial x^2} \tag{2.8}
$$

### 2.5 截面内力

**弯矩的定义与计算**

弯矩$M$定义为轴向应力对中性轴的力矩：

$$
M = \int_A \sigma_{xx} \cdot z \, dA = \int_A \left( -Ez \frac{\partial^2 w}{\partial x^2} \right) z \, dA \tag{2.9}
$$

由于$E$和$\frac{\partial^2 w}{\partial x^2}$与截面位置$z$无关，可以提出积分：

$$
M = -E \frac{\partial^2 w}{\partial x^2} \int_A z^2 \, dA \tag{2.10}
$$

截面惯性矩的定义为：
$$
I = \int_A z^2 \, dA \tag{2.11}
$$

因此：
$$
M = -EI \frac{\partial^2 w}{\partial x^2} \tag{2.12a}
$$

**符号说明**：这里的负号取决于符号约定。若规定使下部纤维受拉的弯矩为正，则：
$$
\boxed{M = EI \frac{\partial^2 w}{\partial x^2}} \tag{2.12b}
$$

这就是著名的**弯矩-曲率关系**。

**剪力的计算**

剪力$V$可以通过弯矩对$x$求导得到（稍后从平衡方程推导）。

### 2.6 运动方程推导（微元平衡法）

取梁上长度为$dx$的微元，分析其受力情况。

**微元受力分析**：

在位置$x$处：
- 左端面：弯矩$M$，剪力$V$
- 右端面：弯矩$M + \frac{\partial M}{\partial x}dx$，剪力$V + \frac{\partial V}{\partial x}dx$
- 分布惯性力：$\rho A \frac{\partial^2 w}{\partial t^2} dx$（向下为正）

**横向力平衡**（$z$方向）：

$$
V - \left( V + \frac{\partial V}{\partial x}dx \right) + \rho A \frac{\partial^2 w}{\partial t^2} dx \cdot (-1) = 0 \tag{2.13}
$$

注意：惯性力$\rho A \frac{\partial^2 w}{\partial t^2}$的方向与加速度方向相反（达朗贝尔原理），这里乘以$(-1)$是因为我们采用$w$向下为正。

简化得到：
$$
\boxed{\frac{\partial V}{\partial x} = \rho A \frac{\partial^2 w}{\partial t^2}} \tag{2.14}
$$

这是**横向平衡方程**。

**力矩平衡**（绕微元中心）：

$$
M - \left( M + \frac{\partial M}{\partial x}dx \right) + V \cdot \frac{dx}{2} + \left( V + \frac{\partial V}{\partial x}dx \right) \cdot \frac{dx}{2} = 0 \tag{2.15}
$$

忽略$(dx)^2$的高阶小量：
$$
-\frac{\partial M}{\partial x}dx + V \cdot dx = 0 \tag{2.16}
$$

即：
$$
\boxed{\frac{\partial M}{\partial x} = V} \tag{2.17}
$$

这是**力矩平衡方程**。注意：在EB理论中忽略了转动惯量，所以右边没有$\rho I \frac{\partial^2 \theta}{\partial t^2}$项。

**消去剪力得到控制方程**

将力矩平衡方程代入横向平衡方程：

$$
\frac{\partial}{\partial x}\left( \frac{\partial M}{\partial x} \right) = \rho A \frac{\partial^2 w}{\partial t^2} \tag{2.18}
$$

即：
$$
\frac{\partial^2 M}{\partial x^2} = \rho A \frac{\partial^2 w}{\partial t^2} \tag{2.19}
$$

将弯矩-曲率关系$M = EI \frac{\partial^2 w}{\partial x^2}$代入：

$$
\frac{\partial^2}{\partial x^2}\left( EI \frac{\partial^2 w}{\partial x^2} \right) = \rho A \frac{\partial^2 w}{\partial t^2} \tag{2.20}
$$

对于等截面梁（$EI$为常数）：

$$
\boxed{EI \frac{\partial^4 w}{\partial x^4} + \rho A \frac{\partial^2 w}{\partial t^2} = 0} \tag{2.21}
$$

这就是**Euler-Bernoulli梁的横向振动控制方程**——一个四阶偏微分方程。

### 2.7 自由振动分析

**分离变量**

设简谐振动解：
$$
w(x,t) = W(x) e^{i\omega t} \tag{2.22}
$$

其中$W(x)$是振型函数（mode shape），$\omega$是圆频率。

代入控制方程：
$$
EI \frac{d^4 W}{dx^4} e^{i\omega t} + \rho A \cdot W(x) \cdot (-\omega^2) e^{i\omega t} = 0 \tag{2.23}
$$

消去$e^{i\omega t}$：
$$
EI \frac{d^4 W}{dx^4} = \rho A \omega^2 W \tag{2.24}
$$

**无量纲化**

引入参数：
$$
\beta^4 = \frac{\rho A \omega^2}{EI} \tag{2.25}
$$

方程变为：
$$
\frac{d^4 W}{dx^4} = \beta^4 W \tag{2.26}
$$

或写成：
$$
\frac{d^4 W}{dx^4} - \beta^4 W = 0 \tag{2.27}
$$

**求解特征方程**

设$W(x) = e^{rx}$，代入得特征方程：
$$
r^4 = \beta^4 \tag{2.28}
$$

解为：
$$
r = \beta, \quad -\beta, \quad i\beta, \quad -i\beta \tag{2.29}
$$

因此通解为：
$$
W(x) = C_1 e^{\beta x} + C_2 e^{-\beta x} + C_3 e^{i\beta x} + C_4 e^{-i\beta x} \tag{2.30}
$$

利用欧拉公式和双曲函数定义，可以改写为更常用的形式：

$$
\boxed{W(x) = C_1 \cosh(\beta x) + C_2 \sinh(\beta x) + C_3 \cos(\beta x) + C_4 \sin(\beta x)} \tag{2.31}
$$

其中：
$$
\begin{aligned}
\cosh(\beta x) &= \frac{e^{\beta x} + e^{-\beta x}}{2} \\
\sinh(\beta x) &= \frac{e^{\beta x} - e^{-\beta x}}{2}
\end{aligned} \tag{2.32}
$$

### 2.8 边界条件

四阶微分方程需要4个边界条件（每端2个）。

**固支端 (Clamped, C)**

物理含义：位移和转角均为零。

$$
\begin{aligned}
w &= 0 \quad &\text{（位移约束）} \\
\frac{\partial w}{\partial x} &= 0 \quad &\text{（转角约束）}
\end{aligned} \tag{2.33}
$$

**简支端 (Pinned/Simply-supported, P)**

物理含义：位移为零，但可以自由转动；弯矩为零。

$$
\begin{aligned}
w &= 0 \quad &\text{（位移约束）} \\
M = EI\frac{\partial^2 w}{\partial x^2} &= 0 \quad &\text{（弯矩为零）}
\end{aligned} \tag{2.34}
$$

**自由端 (Free, F)**

物理含义：无约束，弯矩和剪力均为零。

$$
\begin{aligned}
M = EI\frac{\partial^2 w}{\partial x^2} &= 0 \quad &\text{（弯矩为零）} \\
V = EI\frac{\partial^3 w}{\partial x^3} &= 0 \quad &\text{（剪力为零）}
\end{aligned} \tag{2.35}
$$

**滑动端 (Roller/Sliding, R)**

物理含义：转角为零（截面不能转动），但可以沿横向滑动；剪力为零。

$$
\begin{aligned}
\frac{\partial w}{\partial x} &= 0 \quad &\text{（转角约束）} \\
V = EI\frac{\partial^3 w}{\partial x^3} &= 0 \quad &\text{（剪力为零）}
\end{aligned} \tag{2.36}
$$

### 2.9 简支梁算例（EB理论）

以两端简支梁为例，演示如何求解固有频率。

**边界条件**：
- $x = 0$: $W(0) = 0$, $W''(0) = 0$
- $x = L$: $W(L) = 0$, $W''(L) = 0$

**应用边界条件**：

在$x = 0$处：

$$
W(0) = C_1 \cdot 1 + C_2 \cdot 0 + C_3 \cdot 1 + C_4 \cdot 0 = C_1 + C_3 = 0 \tag{2.37}
$$

$$
W''(0) = C_1 \beta^2 + C_3 (-\beta^2) = \beta^2(C_1 - C_3) = 0 \tag{2.38}
$$

从这两个方程得：$C_1 = C_3 = 0$

在$x = L$处：

$$
W(L) = C_2 \sinh(\beta L) + C_4 \sin(\beta L) = 0 \tag{2.39}
$$

$$
W''(L) = C_2 \beta^2 \sinh(\beta L) - C_4 \beta^2 \sin(\beta L) = 0 \tag{2.40}
$$

从第二个方程：$C_2 \sinh(\beta L) = C_4 \sin(\beta L)$

代入第一个方程：$2C_4 \sin(\beta L) = 0$

由于$C_4 \neq 0$（否则为零解），必须有：

$$
\sin(\beta L) = 0 \tag{2.41}
$$

因此：
$$
\beta_n L = n\pi, \quad n = 1, 2, 3, \ldots \tag{2.42}
$$

**固有频率**：

$$
\beta_n = \frac{n\pi}{L} \tag{2.43}
$$

由$\beta^4 = \frac{\rho A \omega^2}{EI}$：

$$
\omega_n = \beta_n^2 \sqrt{\frac{EI}{\rho A}} = \frac{n^2 \pi^2}{L^2} \sqrt{\frac{EI}{\rho A}} \tag{2.44}
$$

换算为频率（Hz）：

$$
\boxed{f_n = \frac{\omega_n}{2\pi} = \frac{n^2 \pi}{2L^2} \sqrt{\frac{EI}{\rho A}}} \tag{2.45}
$$

**振型函数**：

$$
W_n(x) = C \sin\left( \frac{n\pi x}{L} \right) \tag{2.46}
$$

---

## 3. Timoshenko梁基础

本章介绍Timoshenko梁的理论基础，包括基本假设、运动学关系和本构关系。

### 3.1 理论动机

Euler-Bernoulli理论有两个主要局限：

1. **忽略剪切变形**：实际上，剪力会导致截面发生剪切畸变，使得截面不再垂直于中性轴。
2. **忽略转动惯量**：截面绕其形心转动时具有转动动能，这在高频振动时不可忽略。

Timoshenko在1921年提出了修正理论，将这两个效应纳入考虑。

### 3.2 基本假设

Timoshenko梁理论的假设与EB理论的关键区别：

**假设1：修正的平截面假设**

> 变形前垂直于梁中性轴的平面截面，在变形后仍然保持为平面，但**不再垂直于变形后的中性轴**。

这意味着：
- 截面可以独立于中性轴转动
- 截面转角$\psi$是一个独立变量，不等于挠度斜率$\frac{\partial w}{\partial x}$

**假设2-4**：与EB理论相同（小变形、线弹性、忽略轴向惯性）

**假设5：考虑转动惯量**

> 截面绕其形心转动的动能不可忽略，需要在动能表达式中包含转动惯量项。

### 3.3 运动学关系（位移场）

**独立变量**

在Timoshenko理论中，有两个独立的场变量：
1. $w(x,t)$：中性轴的横向位移
2. $\psi(x,t)$：截面的转角（不等于$\partial w/\partial x$）

**轴向位移**

距中性轴距离为$z$的点，其轴向位移由截面转动引起：
$$
u(x,z,t) = -z \cdot \psi(x,t) \tag{3.1}
$$

注意与EB理论的区别：这里用的是独立的截面转角$\psi$，而非$\frac{\partial w}{\partial x}$。

**剪切角（剪切变形）的定义**

定义剪切角$\gamma$为挠度斜率与截面转角之差：

$$
\boxed{\gamma(x,t) = \frac{\partial w}{\partial x} - \psi} \tag{3.2}
$$

**物理意义**：

- $\frac{\partial w}{\partial x}$：中性轴的斜率（中性轴的切线方向）
- $\psi$：截面的转角（截面法线与竖直方向的夹角）
- $\gamma$：两者之差，即截面法线与中性轴切线之间的夹角

当$\gamma = 0$时，截面垂直于中性轴，退化为EB理论。

### 3.4 应变-位移关系

**弯曲应变（轴向应变）**

$$
\varepsilon_{xx} = \frac{\partial u}{\partial x} = \frac{\partial}{\partial x}(-z\psi) = -z \frac{\partial \psi}{\partial x} \tag{3.3}
$$

与EB理论的区别：EB理论中$\varepsilon_{xx} = -z \frac{\partial^2 w}{\partial x^2}$，而TB理论中是$-z \frac{\partial \psi}{\partial x}$。

**剪切应变**

$$
\gamma_{xz} = \frac{\partial u}{\partial z} + \frac{\partial w}{\partial x} = -\psi + \frac{\partial w}{\partial x} = \frac{\partial w}{\partial x} - \psi = \gamma \tag{3.4}
$$

**重要结论**：在Timoshenko理论中，剪切应变$\gamma_{xz}$不为零，等于我们定义的剪切角$\gamma$。

### 3.5 本构关系（应力-应变关系）

**轴向应力**

$$
\sigma_{xx} = E \varepsilon_{xx} = -Ez \frac{\partial \psi}{\partial x} \tag{3.5}
$$

**剪切应力**

$$
\tau_{xz} = G \gamma_{xz} = G \gamma = G \left( \frac{\partial w}{\partial x} - \psi \right) \tag{3.6}
$$

### 3.6 截面内力

**弯矩**

$$
M = \int_A \sigma_{xx} \cdot z \, dA = \int_A \left( -Ez \frac{\partial \psi}{\partial x} \right) z \, dA = -E \frac{\partial \psi}{\partial x} \int_A z^2 \, dA \tag{3.7}
$$

$$
\boxed{M = EI \frac{\partial \psi}{\partial x}} \tag{3.8}
$$

（这里采用使下部受拉为正弯矩的符号约定）

**剪力**

剪力由剪切应力在截面上的积分得到。然而，实际的剪切应力分布是非均匀的（例如矩形截面为抛物线分布）。为了处理这个问题，引入**剪切修正系数**$\kappa$：

$$
V = \kappa \int_A \tau_{xz} \, dA = \kappa G A \gamma = \kappa G A \left( \frac{\partial w}{\partial x} - \psi \right) \tag{3.9}
$$

$$
\boxed{V = \kappa G A \left( \frac{\partial w}{\partial x} - \psi \right)} \tag{3.10}
$$

### 3.7 剪切修正系数

**$\kappa$的物理意义**

Timoshenko 梁的剪切修正系数用于修正以下两个效应：

1. **剪切应力非均匀分布**：真实的剪切应力分布不是均匀的
2. **翘曲效应**：截面实际上会发生轻微翘曲

常用截面的$\kappa$值（根据Hutchinson 2001）：

| 截面类型 | $\kappa$表达式 | $\nu=0.3$时的 |
|----------|----------------|----------------:|
| 圆形 | $\frac{6(1+\nu)^2}{7+12\nu+4\nu^2}$ | 0.925 |
| 薄壁圆管 | $\frac{1+\nu}{2+\nu}$ | 0.565 |

---

## 4. Timoshenko梁控制方程

本章推导Timoshenko梁的运动控制方程，并将其转化为无量纲形式。

### 4.1 运动方程推导（微元平衡法）

取梁上长度为$dx$的微元进行分析。

**正方向约定说明**：
> 建立微元平衡方程时，采用以下约定：
> - 微元左截面：$V$向上为正，$M$逆时针为正
> - 微元右截面：$V$向下为正，$M$顺时针为正
> - 惯性力采用牛顿第二定律$\sum F = ma$的形式
> - 转角$\psi$逆时针为正

**横向平衡（z方向）**

考虑微元在$z$方向的力平衡：

$$
\sum F_z = 0: \quad V - \left( V + \frac{\partial V}{\partial x}dx \right) = -\rho A \frac{\partial^2 w}{\partial t^2} dx \tag{4.1}
$$

注意：这里我们用的是牛顿第二定律$\sum F = ma$的形式，而非达朗贝尔原理。

简化得：

$$
-\frac{\partial V}{\partial x}dx = -\rho A \frac{\partial^2 w}{\partial t^2} dx \tag{4.2}
$$

$$
\boxed{\frac{\partial V}{\partial x} = \rho A \frac{\partial^2 w}{\partial t^2}} \tag{4.3}
$$

这是**横向运动方程**。

**转动平衡（绕微元中心）**

考虑微元绕其中心的力矩平衡。**关键区别**：在Timoshenko理论中，必须考虑截面的转动惯量。

$$
\sum M_{\text{center}} = J_{center} \cdot \ddot{\theta} \tag{4.4}
$$

微元的转动惯量为：
$$
J_{center} = \rho I \cdot dx \tag{4.5}
$$

力矩平衡方程（绕微元中心，逆时针为正）：

$$
M - \left( M + \frac{\partial M}{\partial x}dx \right) + V \cdot \frac{dx}{2} + \left( V + \frac{\partial V}{\partial x}dx \right) \cdot \frac{dx}{2} = \rho I \cdot dx \cdot \frac{\partial^2 \psi}{\partial t^2} \tag{4.6}
$$

展开并忽略$(dx)^2$的高阶小量：

$$
-\frac{\partial M}{\partial x}dx + V \cdot dx = \rho I \frac{\partial^2 \psi}{\partial t^2} dx \tag{4.7}
$$

整理得：

$$
\boxed{V - \frac{\partial M}{\partial x} = \rho I \frac{\partial^2 \psi}{\partial t^2}} \tag{4.8}
$$

这是**转动运动方程**。注意右边多了转动惯量项$\rho I \frac{\partial^2 \psi}{\partial t^2}$，这是与EB理论的关键区别。

### 4.2 控制方程组

将本构关系代入运动方程：

**代入横向运动方程**

$$
V = \kappa G A \left( \frac{\partial w}{\partial x} - \psi \right) \tag{4.9}
$$

$$
\frac{\partial V}{\partial x} = \kappa G A \left( \frac{\partial^2 w}{\partial x^2} - \frac{\partial \psi}{\partial x} \right) = \rho A \frac{\partial^2 w}{\partial t^2} \tag{4.10}
$$

整理得到**第一个控制方程**：

$$
\boxed{\kappa G A \left( \frac{\partial^2 w}{\partial x^2} - \frac{\partial \psi}{\partial x} \right) = \rho A \frac{\partial^2 w}{\partial t^2}} \tag{4.11}
$$

**代入转动运动方程**

$$
M = EI \frac{\partial \psi}{\partial x} \tag{4.12}
$$

$$
\frac{\partial M}{\partial x} = EI \frac{\partial^2 \psi}{\partial x^2} \tag{4.13}
$$

代入转动运动方程：

$$
\kappa G A \left( \frac{\partial w}{\partial x} - \psi \right) - EI \frac{\partial^2 \psi}{\partial x^2} = \rho I \frac{\partial^2 \psi}{\partial t^2} \tag{4.14}
$$

整理得到**第二个控制方程**：

$$
\boxed{EI \frac{\partial^2 \psi}{\partial x^2} + \kappa G A \left( \frac{\partial w}{\partial x} - \psi \right) = \rho I \frac{\partial^2 \psi}{\partial t^2}} \tag{4.15}
$$

### 4.3 控制方程组总结

Timoshenko梁的控制方程是两个耦合的二阶偏微分方程：

$$
\begin{cases}
\displaystyle \kappa G A \left( \frac{\partial^2 w}{\partial x^2} - \frac{\partial \psi}{\partial x} \right) = \rho A \frac{\partial^2 w}{\partial t^2} \\[12pt]
\displaystyle EI \frac{\partial^2 \psi}{\partial x^2} + \kappa G A \left( \frac{\partial w}{\partial x} - \psi \right) = \rho I \frac{\partial^2 \psi}{\partial t^2}
\end{cases} \tag{4.16}
$$

**与参考文献对照**：上述方程与Khasawneh & Segalman (2019) Eq.(8-9)完全一致。

与EB理论的比较：
- EB理论：1个四阶方程，1个未知函数$w$
- TB理论：2个二阶耦合方程，2个未知函数$w$和$\psi$

### 4.4 无量纲形式（Khasawneh 2019）

为了便于分析和数值计算，引入无量纲变量。

**无量纲变量定义**

设梁长为$L$，定义：

$$
\xi = \frac{x}{L}, \quad \bar{w} = \frac{w}{L}, \quad \tau = \frac{t}{T} \tag{4.17}
$$

其中特征时间$T$定义为：

$$
T = L^2 \sqrt{\frac{\rho A}{EI}} \tag{4.18}
$$

**物理意义**：$T$是EB梁的特征振动周期量级。

**无量纲参数定义**

引入三个关键无量纲参数（Khasawneh 2019, Eq.6-8）：

$$
\boxed{
\begin{aligned}
\alpha &= \frac{A L^2}{I} \quad &\text{（细长比参数）} \\[8pt]
\beta &= \frac{\kappa G A L^2}{EI} \quad &\text{（剪切刚度参数）} \\[8pt]
\gamma &= \frac{\beta}{\alpha} = \frac{\kappa G}{E} \cdot \frac{I}{A} \cdot \frac{A}{I} = \frac{\kappa G}{E} \quad &\text{（材料参数）}
\end{aligned}
} \tag{4.19}
$$

**参数的物理意义**

- **$\alpha$（细长比参数）**：
  - 对于矩形截面：$\alpha = \frac{bh \cdot L^2}{bh^3/12} = \frac{12L^2}{h^2}$
  - $\alpha$越大，梁越细长，剪切效应越弱

- **$\beta$（剪切刚度参数）**：
  - 表示剪切刚度$\kappa GA$与弯曲刚度$EI/L^2$的相对大小
  - $\beta$越大，剪切刚度相对越强

- **$\gamma$（材料参数）**：
  - 只与材料性质有关，与几何无关
  - 对于各向同性材料：$G = \frac{E}{2(1+\nu)}$，故$\gamma = \frac{\kappa}{2(1+\nu)}$
  - 典型值：钢材（$\nu=0.3$）约为0.33

**无量纲控制方程**

对原方程进行无量纲化：

第一个方程除以$\kappa GA/L$：

$$
\frac{\partial^2 \bar{w}}{\partial \xi^2} - \frac{\partial \psi}{\partial \xi} = \frac{\rho A L^2}{\kappa GA} \cdot \frac{1}{T^2} \frac{\partial^2 \bar{w}}{\partial \tau^2} \tag{4.20}
$$

注意到：
$$
\frac{\rho A L^2}{\kappa GA} \cdot \frac{1}{T^2} = \frac{\rho A L^2}{\kappa GA} \cdot \frac{EI}{\rho A L^4} = \frac{EI}{\kappa GA L^2} = \frac{\alpha}{\beta} \tag{4.21}
$$

因此第一个方程变为：

$$
\beta \left( \frac{\partial^2 \bar{w}}{\partial \xi^2} - \frac{\partial \psi}{\partial \xi} \right) = \alpha \frac{\partial^2 \bar{w}}{\partial \tau^2} \tag{4.22}
$$

类似地，第二个方程除以$EI/L^2$并无量纲化：

$$
\frac{\partial^2 \psi}{\partial \xi^2} + \beta \left( \frac{\partial \bar{w}}{\partial \xi} - \psi \right) = \frac{\partial^2 \psi}{\partial \tau^2} \tag{4.23}
$$

**无量纲控制方程组**

$$
\boxed{
\begin{cases}
\displaystyle \beta \left( \frac{\partial^2 \bar{w}}{\partial \xi^2} - \frac{\partial \psi}{\partial \xi} \right) = \alpha \frac{\partial^2 \bar{w}}{\partial \tau^2} \\[12pt]
\displaystyle \frac{\partial^2 \psi}{\partial \xi^2} + \beta \left( \frac{\partial \bar{w}}{\partial \xi} - \psi \right) = \frac{\partial^2 \psi}{\partial \tau^2}
\end{cases}
} \tag{4.24}
$$

---

## 5. Timoshenko梁边界条件与模态分析

本章讨论Timoshenko梁的边界条件类型、自由振动特征方程以及数值稳定性问题。

### 5.1 边界条件

Timoshenko梁的边界条件涉及4个物理量：
1. $w$：横向位移
2. $\psi$：截面转角
3. $M = EI \frac{\partial \psi}{\partial x}$：弯矩
4. $V = \kappa GA \left( \frac{\partial w}{\partial x} - \psi \right)$：剪力

每端需要2个边界条件，具体取决于边界类型。

**固支端 (Clamped, C)**

物理含义：位移和截面转角都被约束。

$$
\boxed{
\begin{aligned}
w &= 0 \\
\psi &= 0
\end{aligned}
} \tag{5.1}
$$

**简支端 (Pinned, P)**

物理含义：位移被约束，但截面可以自由转动；弯矩为零。

$$
\boxed{
\begin{aligned}
w &= 0 \\
M = EI\frac{\partial \psi}{\partial x} &= 0 \quad \Rightarrow \quad \frac{\partial \psi}{\partial x} = 0
\end{aligned}
} \tag{5.2}
$$

**注意**：TB梁的简支条件是$\psi' = 0$，不是$\psi = 0$！这是因为弯矩$M = EI\psi'$。

**自由端 (Free, F)**

物理含义：无约束，弯矩和剪力都为零。

$$
\boxed{
\begin{aligned}
M = EI\frac{\partial \psi}{\partial x} &= 0 \quad \Rightarrow \quad \frac{\partial \psi}{\partial x} = 0 \\
V = \kappa GA\left(\frac{\partial w}{\partial x} - \psi\right) &= 0 \quad \Rightarrow \quad \frac{\partial w}{\partial x} = \psi
\end{aligned}
} \tag{5.3}
$$

**滑动端 (Roller, R)**

物理含义：截面转角被约束（不能转动），但可以沿横向滑动；剪力为零。

$$
\boxed{
\begin{aligned}
\psi &= 0 \\
V = \kappa GA\left(\frac{\partial w}{\partial x} - \psi\right) &= 0 \quad \Rightarrow \quad \frac{\partial w}{\partial x} = \psi = 0
\end{aligned}
} \tag{5.4}
$$

即$\psi = 0$且$\frac{\partial w}{\partial x} = 0$。

### 5.2 边界条件对照表

**边界条件对照表（Khasawneh 2019 Table 3）**

| 边界类型 | 条件1 | 条件2 | 物理意义 |
|----------|-------|-------|----------|
| **Clamped (C)** | $w = 0$ | $\psi = 0$ | 完全固定 |
| **Pinned (P)** | $w = 0$ | $\psi' = 0$ | 铰支，弯矩为零 |
| **Free (F)** | $\psi' = 0$ | $w' - \psi = 0$ | 弯矩、剪力为零 |
| **Roller (R)** | $\psi = 0$ | $w' - \psi = 0$ | 截面不转，剪力为零 |

### 5.3 自由振动特征方程

**设简谐解**

设时间依赖为简谐形式：

$$
\bar{w}(\xi,\tau) = \hat{w}(\xi) e^{i\sqrt{\lambda}\tau}, \quad \psi(\xi,\tau) = \hat{\psi}(\xi) e^{i\sqrt{\lambda}\tau} \tag{5.5}
$$

其中$\lambda$是无量纲特征值。

**与物理频率的关系**：

$$
\lambda = \omega^2 T^2 = \omega^2 L^4 \frac{\rho A}{EI} \tag{5.6}
$$

因此：
$$
\omega = \frac{\sqrt{\lambda}}{L^2} \sqrt{\frac{EI}{\rho A}} \tag{5.7}
$$

**代入控制方程**

代入后消去$e^{i\sqrt{\lambda}\tau}$：

$$
\begin{cases}
\displaystyle \beta \left( \frac{d^2 \hat{w}}{d\xi^2} - \frac{d\hat{\psi}}{d\xi} \right) = -\alpha \lambda \hat{w} \\[12pt]
\displaystyle \frac{d^2 \hat{\psi}}{d\xi^2} + \beta \left( \frac{d\hat{w}}{d\xi} - \hat{\psi} \right) = -\lambda \hat{\psi}
\end{cases} \tag{5.8}
$$

### 5.4 特征根分析

设解的形式为$\hat{w} = W_0 e^{r\xi}$，$\hat{\psi} = \Psi_0 e^{r\xi}$，代入方程组可以得到关于$r$的特征方程。

特征根的性质取决于$\lambda$与$\alpha$的关系：

**情况1：$\lambda < \alpha$（低频模态）**

此时特征方程有四个根：两个实根$\pm p$和两个纯虚根$\pm iq$。

通解形式：
$$
\hat{w}(\xi) = A_1 \cosh(p\xi) + A_2 \sinh(p\xi) + A_3 \cos(q\xi) + A_4 \sin(q\xi) \tag{5.9}
$$

其中：
$$
\begin{aligned}
p^2 &= \frac{1}{2}\left[ \lambda\left(1 + \frac{\alpha}{\beta}\right) + \sqrt{\lambda^2\left(1 - \frac{\alpha}{\beta}\right)^2 + 4\alpha\beta} \right] - \beta \\
q^2 &= -\frac{1}{2}\left[ \lambda\left(1 + \frac{\alpha}{\beta}\right) - \sqrt{\lambda^2\left(1 - \frac{\alpha}{\beta}\right)^2 + 4\alpha\beta} \right] + \beta
\end{aligned} \tag{5.10}
$$

**情况2：$\lambda = \alpha$（临界频率）**

此时特征方程退化，需要特殊处理。这个频率称为**第二频谱截止频率**。

**情况3：$\lambda > \alpha$（高频模态）**

此时所有四个特征根都可能是复数，模态解的形式更加复杂。

详细推导见Khasawneh & Segalman (2019) Table 4-8。

### 5.5 数值稳定性问题与解决方案

**问题描述**

在求解高阶模态时，传统的解析解形式会遇到数值不稳定问题：

1. **双曲函数溢出**：$\cosh(\beta L)$和$\sinh(\beta L)$在$\beta L$较大时会溢出
2. **有效数字损失**：大数相减导致精度丢失

**Khasawneh & Segalman (2019)的解决方案**

引入稳定化函数：

$$
f(\mu) = \frac{\sinh(\mu)}{\mu}, \quad g(\mu) = \frac{\cosh(\mu) - 1}{\mu^2}, \quad \psi(\mu) = \frac{\sinh(\mu) - \mu}{\mu^3} \tag{5.11}
$$

这些函数在$\mu \to 0$时有良好的极限行为，避免了0/0型不定式。

详细的稳定化特征方程见Khasawneh & Segalman (2019) Table 8和Table 12。

---

## 6. 两种理论的详细对比

### 6.1 基本假设对比

| 项目 | Euler-Bernoulli | Timoshenko |
|------|-----------------|------------|
| **截面变形** | 保持垂直于中性轴 | 保持平面但可倾斜 |
| **剪切变形** | $\gamma = 0$（忽略） | $\gamma = w' - \psi \neq 0$ |
| **转动惯量** | 忽略 | 包含$\rho I \ddot{\psi}$ |
| **独立变量** | 1个（$w$） | 2个（$w$, $\psi$） |
| **本构关系** | $M = EI w''$ | $M = EI\psi'$, $V = \kappa GA(w'-\psi)$ |

### 6.2 控制方程对比

| 项目 | Euler-Bernoulli | Timoshenko |
|------|-----------------|------------|
| **方程数** | 1 | 2（耦合） |
| **方程阶数** | 4阶 | 各2阶 |
| **方程形式** | $EIw'''' + \rho A\ddot{w} = 0$ | 见第4章方程组 |
| **边界条件数** | 4（每端2个） | 4（每端2个） |

### 6.3 频率关系

对于给定的梁，Timoshenko频率总是**低于或等于**Euler-Bernoulli频率：

$$
\omega_n^{TB} \leq \omega_n^{EB} \tag{6.1}
$$

**物理解释**：
- TB理论考虑了额外的变形模式（剪切变形）
- 转动惯量增加了系统的惯性

**频率比的近似公式**（Rayleigh修正）：

对于简支梁第$n$阶模态：

$$
\frac{\omega_n^{TB}}{\omega_n^{EB}} \approx \frac{1}{\sqrt{1 + \frac{n^2\pi^2}{L^2}\left(\frac{I}{A} + \frac{EI}{\kappa GA}\right)}} \tag{6.2}
$$

记$r^2 = I/A$（回转半径的平方），$s^2 = EI/(\kappa GA)$（剪切变形参数），则：

$$
\frac{\omega_n^{TB}}{\omega_n^{EB}} \approx \frac{1}{\sqrt{1 + \frac{n^2\pi^2(r^2 + s^2)}{L^2}}} \tag{6.3}
$$

### 6.4 适用范围

**Euler-Bernoulli理论适用条件**：

1. 细长比：$L/h > 10$（或更保守地$L/h > 20$）
2. 只关心低阶模态（通常前3-5阶）
3. 各向同性材料

**必须使用Timoshenko理论的情况**：

1. 短粗梁：$L/h < 10$
2. 高阶模态分析
3. 复合材料梁（$G/E$比值较小）
4. 夹层结构
5. 高精度要求的场合

### 6.5 误差估计

相对频率误差：

$$
\epsilon = \frac{\omega_{EB} - \omega_{TB}}{\omega_{TB}} \times 100\% \tag{6.4}
$$

对于矩形截面梁（$\nu = 0.3$）：

| $L/h$ | 第1阶误差 | 第5阶误差 | 第10阶误差 |
|-------|-----------|-----------|------------|
| 20 | 0.2% | 5% | 19% |
| 10 | 0.8% | 17% | 53% |
| 5 | 3% | 45% | 超过100% |

**结论**：即使是细长梁，高阶模态也需要使用Timoshenko理论。

---

## 7. Hamilton原理推导（能量方法）

除了微元平衡法，控制方程也可以通过Hamilton原理（变分法）推导。这种方法对于有限元推导尤为重要。

### 7.1 Hamilton原理陈述

对于保守系统：

$$
\delta \int_{t_1}^{t_2} (T - U) \, dt = 0 \tag{7.1}
$$

其中$T$是动能，$U$是势能（应变能）。

**变分原理的含义**：真实运动路径使作用量泛函取驻值。

### 7.2 动能表达式

Timoshenko梁的动能包含两部分：

**平动动能**（横向运动）：
$$
T_w = \frac{1}{2} \int_0^L \rho A \left( \frac{\partial w}{\partial t} \right)^2 dx \tag{7.2}
$$

**转动动能**（截面转动）：
$$
T_\psi = \frac{1}{2} \int_0^L \rho I \left( \frac{\partial \psi}{\partial t} \right)^2 dx \tag{7.3}
$$

**总动能**：
$$
\boxed{T = \frac{1}{2} \int_0^L \left[ \rho A \left( \frac{\partial w}{\partial t} \right)^2 + \rho I \left( \frac{\partial \psi}{\partial t} \right)^2 \right] dx} \tag{7.4}
$$

### 7.3 势能（应变能）表达式

应变能包含弯曲应变能和剪切应变能：

**弯曲应变能**：
$$
U_b = \frac{1}{2} \int_0^L EI \left( \frac{\partial \psi}{\partial x} \right)^2 dx \tag{7.5}
$$

**剪切应变能**：
$$
U_s = \frac{1}{2} \int_0^L \kappa GA \left( \frac{\partial w}{\partial x} - \psi \right)^2 dx \tag{7.6}
$$

**总势能**：
$$
\boxed{U = \frac{1}{2} \int_0^L \left[ EI \left( \frac{\partial \psi}{\partial x} \right)^2 + \kappa GA \left( \frac{\partial w}{\partial x} - \psi \right)^2 \right] dx} \tag{7.7}
$$

### 7.4 变分推导详解

拉格朗日函数：
$$
L = T - U \tag{7.8}
$$

Hamilton原理要求：
$$
\delta \int_{t_1}^{t_2} L \, dt = \delta \int_{t_1}^{t_2} (T - U) \, dt = 0 \tag{7.9}
$$

#### 7.4.1 对$w$的变分

**步骤1：动能变分**

$$
\delta T_w = \int_0^L \rho A \dot{w} \, \delta\dot{w} \, dx \tag{7.10}
$$

对时间积分，应用分部积分：

$$
\int_{t_1}^{t_2} \delta T_w \, dt = \int_{t_1}^{t_2} \int_0^L \rho A \dot{w} \, \delta\dot{w} \, dx \, dt \tag{7.11}
$$

交换积分顺序，对时间分部积分：

$$
\int_{t_1}^{t_2} \rho A \dot{w} \, \delta\dot{w} \, dt = \left[ \rho A \dot{w} \, \delta w \right]_{t_1}^{t_2} - \int_{t_1}^{t_2} \rho A \ddot{w} \, \delta w \, dt \tag{7.12}
$$

边界条件要求$\delta w(t_1) = \delta w(t_2) = 0$，因此边界项为零：

$$
\int_{t_1}^{t_2} \delta T_w \, dt = -\int_{t_1}^{t_2} \int_0^L \rho A \ddot{w} \, \delta w \, dx \, dt \tag{7.13}
$$

**步骤2：剪切应变能变分**

$$
\delta U_s = \int_0^L \kappa GA \left( \frac{\partial w}{\partial x} - \psi \right) \left( \frac{\partial \delta w}{\partial x} \right) dx \tag{7.14}
$$

对空间分部积分：

$$
\int_0^L \kappa GA (w' - \psi) \, \delta w' \, dx = \left[ \kappa GA (w' - \psi) \, \delta w \right]_0^L - \int_0^L \kappa GA (w'' - \psi') \, \delta w \, dx \tag{7.15}
$$

**步骤3：组合变分结果**

将以上结果代入Hamilton原理：

$$
\int_{t_1}^{t_2} \left\{ -\int_0^L \rho A \ddot{w} \, \delta w \, dx - \left[ \kappa GA (w' - \psi) \, \delta w \right]_0^L + \int_0^L \kappa GA (w'' - \psi') \, \delta w \, dx \right\} dt = 0 \tag{7.16}
$$

由$\delta w$的任意性（在区域内部），得到**第一个控制方程**：

$$
\boxed{\kappa GA (w'' - \psi') = \rho A \ddot{w}} \tag{7.17}
$$

边界项给出**自然边界条件**（剪力条件）：

$$
V = \kappa GA (w' - \psi) = 0 \quad \text{或} \quad w \text{ 给定} \tag{7.18}
$$

#### 7.4.2 对$\psi$的变分

**步骤1：转动动能变分**

$$
\delta T_\psi = \int_0^L \rho I \dot{\psi} \, \delta\dot{\psi} \, dx \tag{7.19}
$$

对时间分部积分（类似上面的过程）：

$$
\int_{t_1}^{t_2} \delta T_\psi \, dt = -\int_{t_1}^{t_2} \int_0^L \rho I \ddot{\psi} \, \delta\psi \, dx \, dt \tag{7.20}
$$

**步骤2：弯曲应变能变分**

$$
\delta U_b = \int_0^L EI \psi' \, \delta\psi' \, dx \tag{7.21}
$$

对空间分部积分：

$$
\int_0^L EI \psi' \, \delta\psi' \, dx = \left[ EI \psi' \, \delta\psi \right]_0^L - \int_0^L EI \psi'' \, \delta\psi \, dx \tag{7.22}
$$

**步骤3：剪切应变能对$\psi$的变分**

$$
\delta_\psi U_s = \int_0^L \kappa GA (w' - \psi) \cdot (-\delta\psi) \, dx = -\int_0^L \kappa GA (w' - \psi) \, \delta\psi \, dx \tag{7.23}
$$

**步骤4：组合变分结果**

$$
\int_{t_1}^{t_2} \left\{ -\int_0^L \rho I \ddot{\psi} \, \delta\psi \, dx - \left[ EI \psi' \, \delta\psi \right]_0^L + \int_0^L EI \psi'' \, \delta\psi \, dx + \int_0^L \kappa GA (w' - \psi) \, \delta\psi \, dx \right\} dt = 0 \tag{7.24}
$$

由$\delta\psi$的任意性，得到**第二个控制方程**：

$$
\boxed{EI \psi'' + \kappa GA (w' - \psi) = \rho I \ddot{\psi}} \tag{7.25}
$$

边界项给出**自然边界条件**（弯矩条件）：

$$
M = EI \psi' = 0 \quad \text{或} \quad \psi \text{ 给定} \tag{7.26}
$$

### 7.5 变分推导总结

通过Hamilton原理，我们得到了与微元平衡法完全相同的控制方程：

$$
\begin{cases}
\displaystyle \kappa GA (w'' - \psi') = \rho A \ddot{w} \\[8pt]
\displaystyle EI \psi'' + \kappa GA (w' - \psi) = \rho I \ddot{\psi}
\end{cases} \tag{7.27}
$$

**Hamilton原理的优势**：
1. 能量表达式可以直接用于有限元离散
2. 边界条件自然地从变分过程中导出
3. 为更复杂的梁理论（如考虑预应力、复合材料等）提供了统一的推导框架

---

## 8. 参考文献

1. **Khasawneh, F. A., & Segalman, D. J. (2019)**. Exact and numerically stable expressions for Euler-Bernoulli and Timoshenko beam modes. *International Journal of Mechanical Sciences*, 166, 105234.
   - DOI: https://doi.org/10.1016/j.ijmecsci.2019.105234
   - **主要参考内容**：无量纲参数定义、边界条件表、数值稳定解析解

2. **Friedman, Z., & Kosmatka, J. B. (1993)**. An improved two-node Timoshenko beam finite element. *Computers & Structures*, 47(3), 473-481.
   - DOI: https://doi.org/10.1016/0045-7949(93)90243-7
   - **主要参考内容**：Hamilton原理推导、能量表达式

3. **Hutchinson, J. R. (2001)**. Shear coefficients for Timoshenko beam theory. *Journal of Applied Mechanics*, 68(1), 87-92.
   - DOI: https://doi.org/10.1115/1.1349417
   - **主要参考内容**：剪切修正系数

4. **Timoshenko, S. P. (1921)**. On the correction for shear of the differential equation for transverse vibrations of prismatic bars. *Philosophical Magazine*, 41(245), 744-746.
   - **原始文献**：Timoshenko梁理论的首次提出

5. **Han, S. M., Benaroya, H., & Wei, T. (1999)**. Dynamics of transversely vibrating beams using four engineering theories. *Journal of Sound and Vibration*, 225(5), 935-988.
   - **综述参考**：各种梁理论的详细对比

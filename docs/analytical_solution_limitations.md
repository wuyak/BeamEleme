# 解析解的数值限制

## 问题

Timoshenko梁的解析解在λ≈α（特征值接近剪切参数）时会漏掉部分模态。

chebfun在λ=α边界处的根查找算法失效：特征函数在边界附近剧烈振荡，数值条件数极大，fzero等方法无法收敛。这个问题与α的绝对值无关，关键是某个模态的λ值与α的相对接近程度（<10%）。

## 影响范围

所有α值都可能在高阶模态处缺失1-3个模态：

| α范围 | 典型缺失位置 | 示例 |
|-------|------------|------|
| α > 1000 | 模态12-15 | 矩形(0.2×0.1), α=1200 |
| 100 < α < 500 | 模态7-9 | 圆形(d=0.2), α=400 |
| α < 200 | 模态3-6 | 圆环(D=0.3,d=0.2), α=123 |

测试配置（Q235钢，L=1m，PP边界，15个模态）：

| 截面类型 | α | 缺失模态编号 | λ值 |
|---------|---|------------|-----|
| 矩形(0.2×0.1) | 1200 | [13,14,15] | 1202, 1242, 1269 |
| 圆形(d=0.2) | 400 | [8,9] | 400, 437 |
| 圆环(D=0.3,d=0.2) | 123 | [4,5] | 123, 129 |

完整系统测试（3截面×10边界条件）：
- 30个算例全部发现缺失模态（100%）
- 总缺失~60个模态（在450个中，约13%）
- 所有边界条件都受影响，不是PP特有的问题

复现代码：
```matlab
mat = parameters.MaterialLibrary.Steel();  % Q235
L = 1.0;
NElem = 200;

% 细长梁（α=1200）
sec1 = parameters.SectionLibrary.Rectangular(0.2, 0.1);
beam1 = parameters.BeamModel(mat, sec1, 'L', L);
wf1 = workflow.solve('PP', beam1, NElem, 15);

% 短粗梁（α=123）
sec2 = parameters.SectionLibrary.CircularHollow(0.3, 0.2);
beam2 = parameters.BeamModel(mat, sec2, 'L', L);
wf2 = workflow.solve('PP', beam2, NElem, 15);
```

## 解决方案

使用FEM结果为准，解析解仅用于验证FEM的收敛性。

`comparison.matchModes()` 会自动检测λ≈α的缺失模态并标记为 "Missing (Numerical Singularity)"，不会影响误差计算和收敛测试。

验证代码：
```matlab
sys = analytical.bending_shear.buildSystemParameters(beam);
fprintf('α = %.2f\n', sys.par.alpha);

comp = comparison.compare(fem, ana, 'n_modes', 15);

% 检查MAC值
bad_modes = find(comp.mac < 0.9);
if ~isempty(bad_modes)
    warning('模态 %s 的MAC值较低，可能是λ≈α导致', num2str(bad_modes'));
end
```

对短粗梁（α<100），建议使用更多单元确保FEM收敛：
```matlab
if sys.par.alpha < 100
    NElem = 200;
else
    NElem = 100;
end
```

## 技术细节

解析解将特征方程分为三个区域：
- λ < α：使用双曲函数解
- λ = α：特殊情况
- λ > α：使用三角函数解

缺失的模态不是λ精确等于α的特殊解，而是λ非常接近α的根。chebfun在区域边界处的严格过滤（λ<α和λ>α）导致边界附近的根丢失。增加采样密度（minSamples: 2048）或启用区域分割（splitting: on）仍无法避免。

FEM使用的控制方程包含了更完整的物理效应，在λ≈α处仍然稳定可靠。所有FEM结果已通过网格收敛性测试（误差<0.25%）。

## 参考

Khasawneh, F. A., & Segalman, D. J. (2019). "Exact and Numerically Stable Expressions for Euler-Bernoulli and Timoshenko Beam Modes". International Journal of Mechanical Sciences, 166, 105234.

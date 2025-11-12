# References

本项目的理论推导和实现基于以下文献和开源代码。所有文献均可通过 DOI 或 ISBN 在学术数据库中查询获取。

---

## 期刊论文

### Friedman1993
> Friedman, Z., & Kosmatka, J. B. (1993). An improved two-node Timoshenko beam finite element. *Computers & Structures*, 47(3), 473-481.

- **DOI**: https://doi.org/10.1016/0045-7949(93)90243-7
- **引用内容**: 
  - 避免剪切锁定的 Timoshenko 梁精确刚度矩阵（Appendix, Eq. A1）
  - 一致质量矩阵推导（Appendix, Eq. A4–A6）

### Khasawneh2019
> Khasawneh, F. A., & Segalman, D. J. (2019). Exact and numerically stable expressions for Euler-Bernoulli and Timoshenko beam modes. *International Journal of Mechanical Sciences*, 166, 105234.

- **DOI**: https://doi.org/10.1016/j.ijmecsci.2019.105234
- **引用内容**: 
  - Timoshenko 梁边界条件的约束设定（Table 3）
  - Timoshenko 梁弯曲-剪切耦合模态的传统解析解（Eq. 15-19）
  - Timoshenko 梁弯曲-剪切耦合模态的改进稳定解析解（Table 8, Table 12）

### Hutchinson2001
> Hutchinson, J. R. (2001). Shear coefficients for Timoshenko beam theory. *Journal of Applied Mechanics*, 68(1), 87-92.

- **DOI**: https://doi.org/10.1115/1.1349417
- **引用内容**: 
  - Timoshenko 梁矩形截面剪切修正系数（Eq. 46–47）
  - Timoshenko 梁圆形截面剪切修正系数（Eq. 43）
  - Timoshenko 梁环形截面剪切修正系数（Eq. 44, 53）

---

## 专著

### Roark's Formulas
> Young, W. C., Budynas, R. G., & Sadegh, A. M. (2012). *Roark's Formulas for Stress and Strain* (8th ed.). McGraw-Hill.

- **ISBN**: 978-0-07-174248-1
- **引用内容**: 
  - 截面几何性质（Appendix A, Table A.1）：面积、惯性矩
  - 截面抗扭惯性矩 (Chapter 10, Table 10.1)

### Blevins1979
> Blevins, R. D. (1979). *Formulas for Natural Frequency and Mode Shape*. Van Nostrand Reinhold.

- **ISBN**: 0-442-20710-7
- **引用内容**: 
  - 梁的轴向振动频率（Table 8-16）
  - 梁的扭转振动频率（Table 8-19）

---

## 开源代码

### Khasawneh2019Code
> Khasawneh, F. A., & Segalman, D. J. (2018). Code for: Exact and Numerically Stable Expressions for Euler-Bernoulli and Timoshenko Beam Modes [Dataset]. Mendeley Data, V1.

- **获取**: https://doi.org/10.17632/r275tx2yp8.1
- **许可证**: CC BY 4.0
- **引用内容**：Timoshenko梁弯曲-剪切模态分析的解析解实现

**对代码的引用和修正**：详见 [`docs/code_corrections.md`](../docs/code_corrections.md)

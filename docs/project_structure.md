# 项目目录结构

## 完整目录树

```
BeamElem/
├── example_basic.m                     # 基本使用示例
├── chebfun/                            # chebfun 开源计算库
├── .cache/                             # 计算缓存
├── src/                                # 核心源代码
│   ├── +workflow/                      # 求解流程控制
│   │   ├── solve.m
│   │   └── tests/
│   │       ├── test_workflow.m
│   │       └── convergence/            # 收敛率测试
│   │           ├── convergence.m
│   │           └── output/
│   ├── +parameters/                    # 参数管理
│   │   ├── BeamModel.m
│   │   ├── MaterialLibrary.m
│   │   ├── SectionLibrary.m
│   │   └── SectionCalculator.m
│   ├── +solvers/                       # 有限元求解器
│   │   ├── solve.m
│   │   ├── assemble.m
│   │   ├── eigensolve.m
│   │   ├── filterEigen.m
│   │   ├── buildResult.m
│   │   ├── getDOFPerNode.m
│   │   ├── getRigidModeCount.m
│   │   └── handleEigsFailure.m
│   ├── +matrix/                        # 单元矩阵生成
│   │   ├── +stiffness/                 # 刚度矩阵
│   │   │   ├── +axial/
│   │   │   │   └── element.m
│   │   │   ├── +torsion/
│   │   │   │   └── element.m
│   │   │   └── +bending_shear/
│   │   │       └── element.m
│   │   └── +mass/                      # 质量矩阵
│   │       ├── +axial/
│   │       │   └── element.m
│   │       ├── +torsion/
│   │       │   └── element.m
│   │       └── +bending_shear/
│   │           └── element.m
│   ├── +analytical/                    # 解析解模块
│   │   ├── +axial/
│   │   │   ├── solve.m
│   │   │   ├── getFrequencies.m
│   │   │   └── getModes.m
│   │   ├── +torsion/
│   │   │   ├── solve.m
│   │   │   ├── getFrequencies.m
│   │   │   └── getModes.m
│   │   └── +bending_shear/
│   │       ├── solve.m
│   │       ├── getEigenvalues.m
│   │       ├── getModes.m
│   │       └── buildSystemParameters.m
│   ├── +comparison/                    # FEM与解析解对比
│   │   ├── compare.m
│   │   ├── matchModes.m
│   │   └── computeMAC.m
│   ├── +boundary_conditions/           # 边界条件管理
│   │   └── BoundaryConditionManager.m
│   └── +visualization/                 # 可视化工具
│       ├── plotMode.m
│       ├── loadCache.m
│       ├── saveFigure.m
│       └── tests/
│           ├── test_api.m
│           ├── test_complete.m
│           ├── demo.m
│           └── output/
├── reference/                          # 理论基础文献
│   └── REFERENCES.md
└── docs/                               # 项目文档
    ├── project_structure.md            # 本文档
    ├── code_corrections.md
    ├── analytical_solution_limitations.md
    ├── timoshenko_beam_theory.md       # Timoshenko 梁理论推导
    ├── axial_torsion_modes.md          # 轴向与扭转振动模态
    ├── bending_shear_matrices.md       # 刚度矩阵与质量矩阵推导
    └── bending_shear_analytical.md     # 解析解推导
```

## 快速查找

| 需求 | 修改位置 |
|------|---------|
| 添加新材料 | `+parameters/MaterialLibrary.m` |
| 添加新截面类型 | `+parameters/SectionLibrary.m` 和 `SectionCalculator.m` |
| 修改单元矩阵 | `+matrix/+stiffness/` 或 `+matrix/+mass/` |
| 添加新边界条件 | `+boundary_conditions/BoundaryConditionManager.m` |
| 修改解析解算法 | `+analytical/+<振动类型>/` |
| 调整验证逻辑 | `+comparison/` |
| 修改可视化样式 | `+visualization/plotMode.m` |

## 扩展指南

### 简单扩展

#### 1. 添加新材料

**位置**：`+parameters/MaterialLibrary.m`

在 `MaterialLibrary` 类中添加静态方法：

```matlab
function mat = NewMaterial()
    mat = struct(...
        'E',   210e9,  % 弹性模量 (Pa)
        'G',   80e9,   % 剪切模量 (Pa)
        'rho', 7850,   % 密度 (kg/m³)
        'nu',  0.3 );  % 泊松比
end
```

使用：
```matlab
mat = parameters.MaterialLibrary.NewMaterial();
```

#### 2. 添加新截面类型

**位置**：
- `+parameters/SectionLibrary.m` — 添加截面创建方法
- `+parameters/SectionCalculator.m` — 添加截面属性计算逻辑

**Step 1**：在 `SectionLibrary.m` 中添加截面创建方法

```matlab
function sec = NewSection(param1, param2)
    sec = struct(...
        'type', 'NewSection', ...
        'param1', param1, ...
        'param2', param2 );
end
```

**Step 2**：在 `SectionCalculator.m` 的 `calculate()` 方法中添加计算逻辑

```matlab
case 'NewSection'
    A = ...;    % 截面面积
    I = ...;    % 惯性矩
    J = ...;    % 扭转常数
    kappa = ...; % 剪切修正系数
```

### 复杂扩展

#### 3. 扩展到弯曲-扭转耦合

**涉及模块**：`+matrix`、`+solvers`、`+analytical`、`+workflow`

**Step 1：单元矩阵 (+matrix/)**

在 `+matrix`  中新建 `+coupling` ，根据耦合算法调用调整已有的刚度矩阵和质量矩阵。

**Step 2：求解器 (+solvers/)**

修改以下文件以支持耦合模型：
- `assemble.m` — 支持新的耦合矩阵装配
- `getDOFPerNode.m` — 返回正确的自由度数（耦合模型可能是 2 或 4 个自由度）

**Step 3：解析解 (+analytical/)**

如果耦合振动存在闭式解：
- 在 `+analytical/` 下新建 `+bending_torsion/` 子包
- 实现 `solve.m`、`getFrequencies.m`、`getModes.m`

如果不存在闭式解，可跳过此步，仅使用 FEM 求解。

**Step 4：工作流 (+workflow/)**

修改 `workflow.solve()` 添加新的振动类型分支：

**Step 5：可视化 (+visualization/)**

扩展 `plotMode.m` 支持耦合模态的多变量显示（如同时绘制弯曲和扭转位移）。

#### 4. 与其他模型耦合

如需将 Timoshenko 梁与其他结构（如板壳、实体单元）耦合：

**关键步骤**：

1. **接口设计** — 定义梁单元与其他单元的节点共享规则
   - 确定自由度匹配方式（如梁的转角自由度如何映射到实体单元）

2. **矩阵装配** — 修改 `+solvers/assemble.m` 支持混合单元类型
   - 实现多类型单元的刚度矩阵和质量矩阵合成

3. **边界条件** — 扩展 `+boundary_conditions/BoundaryConditionManager.m`
   - 处理复杂约束（如梁端与板边界的连接约束）

4. **验证流程** — 可能需要禁用解析解对比
   - 因为耦合模型通常无闭式解
   - 可采用其他验证方式（如与商业软件对比）

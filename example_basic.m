%% BeamElem 基本使用示例
% 本脚本演示如何定义材料、截面、求解并可视化结果

clear; clc;

%% ============================================================
%  第一步：定义材料参数
%  材料定义位置：src/+parameters/MaterialLibrary.m
%
%  目前内置材料：
%    - Steel()    : Q235 结构钢 (E=200 GPa, ρ=7850 kg/m³, ν=0.3)
%    - Aluminum() : 铝合金 (E=70 GPa, ρ=2700 kg/m³, ν=0.33)
%    - Titanium() : 钛合金 (E=110 GPa, ρ=4500 kg/m³, ν=0.34)
%
%  自定义材料示例：
%    mat = parameters.MaterialProperties(200e9, 7850, 0.3);  % E, ρ, ν
% ============================================================

mat = parameters.MaterialLibrary.Steel();


%% ============================================================
%  第二步：定义截面参数
%  截面定义位置：src/+parameters/SectionLibrary.m
%
%  目前支持的截面类型：
%    1. 矩形截面：Rectangular(width, height)
%       - 要求：width ≥ height（宽度不小于高度）
%       - 示例：SectionLibrary.Rectangular(0.2, 0.1)  % 宽20cm×高10cm
%
%    2. 圆形截面：Circular(diameter)
%       - 示例：SectionLibrary.Circular(0.15)  % 直径15cm
%
%    3. 环形截面：Ring(outer_diameter, inner_diameter)
%       - 要求：outer_diameter > inner_diameter
%       - 示例：SectionLibrary.Ring(0.2, 0.18)  % 外径20cm, 内径18cm
%
%  如需添加新截面类型：
%    1. 在 SectionLibrary.m 中添加静态方法
%    2. 在 SectionCalculator.m 中添加计算公式（A, I, J, κ）
%       - A: 截面积 (参考 Roark's formulas)
%       - I: 惯性矩 (参考 Roark's formulas)
%       - J: 抗扭惯性矩 (参考 Roark's formulas)
%       - κ: 剪切修正系数 (参考 Hutchinson 2000)
% ============================================================

sec = parameters.SectionLibrary.Rectangular(0.2, 0.1);  % 矩形截面：宽20cm×高10cm


%% ============================================================
%  第三步：创建梁模型
%  参数说明：
%    - mat: 材料对象
%    - sec: 截面对象
%    - 'L': 梁长度（单位：米）
% ============================================================

beam = parameters.BeamModel(mat, sec, 'L', 1.0);  % 创建长度为1米的梁


%% ============================================================
%  第四步：求解
%  参数说明：
%    - 边界条件: 两个字符的字符串（第一个字符是左端，第二个是右端）
%      * P (Pinned)  : 铰支 - 约束位移，不约束转角
%      * C (Clamped) : 固支 - 约束位移和转角
%      * R (Roller)  : 滚支 - 约束横向位移
%      * F (Free)    : 自由 - 无约束
%    - beam: 梁模型对象
%    - NElem: 有限元单元数（建议 ≥ 100）
%    - n_modes: 求解的模态数
%
%  支持的边界条件组合（10种）：
%    'PP', 'CF', 'CC', 'FF', 'CP', 'PR', 'CR', 'RR', 'RF', 'PF'
%
%  求解过程：
%    1. 装配有限元刚度和质量矩阵
%    2. 求解特征值问题（频率和模态）
%    3. 计算解析解（用于验证）
%    4. 对比 FEM 和解析解（频率误差、MAC值）
%
%  缓存规则：
%    - 缓存位置：.cache/ 目录
%    - 命名规则：<hash>_<BC>_N<NElem>_M<n_modes>/
%      * hash: 梁模型的MD5哈希值（材料+截面+长度）
%      * BC: 边界条件
%      * NElem: 单元数
%      * n_modes: 模态数
%    - 示例：.cache/4BB0F885_PP_N100_M5/
%    - 如果缓存存在，会自动加载，避免重复计算
% ============================================================

fprintf('开始求解...\n');
wf_result = workflow.solve('PP', beam, 100, 10);  % 简支梁，100单元，10模态


%% ============================================================
%  第五步：查看结果
%  wf_result 包含三个字段：
%    - axial: 轴向振动结果
%    - torsion: 扭转振动结果
%    - bending_shear: 弯曲-剪切振动结果
%
%  每个结果包含：
%    - fem: 有限元结果
%    - ana: 解析解结果
%    - comp: 对比结果（频率误差、MAC值）
% ============================================================

fprintf('\n======== 结果统计 ========\n');

% 轴向振动
comp_axial = wf_result.axial.comp;
fprintf('【轴向振动】\n');
fprintf('  最大频率误差: %.4f%%\n', max(comp_axial.freq_error));
fprintf('  平均MAC值: %.6f\n', mean(comp_axial.mac));

% 扭转振动
comp_torsion = wf_result.torsion.comp;
fprintf('【扭转振动】\n');
fprintf('  最大频率误差: %.4f%%\n', max(comp_torsion.freq_error));
fprintf('  平均MAC值: %.6f\n', mean(comp_torsion.mac));

% 弯曲-剪切振动
comp_bending = wf_result.bending_shear.comp;
fprintf('【弯曲-剪切振动】\n');
fprintf('  最大频率误差: %.4f%%\n', max(comp_bending.freq_error));
fprintf('  平均MAC值: %.6f\n', mean(comp_bending.mac));


%% ============================================================
%  第六步：可视化（可选）
%  使用 visualization 模块绘制模态振型
%
%  visualization.plotMode 的正确用法：
%    h = visualization.plotMode(BC, beam, NElem, mode_idx)
%    h = visualization.plotMode(BC, beam, NElem, mode_idx, 'matrix_method', 'bending_shear')
%
%  参数说明：
%    - BC: 边界条件字符串 ('PP', 'CF', 等)
%    - beam: BeamModel 对象
%    - NElem: 单元数
%    - mode_idx: 模态索引（标量或向量，如 1 或 [1 2 3]）
%    - 'matrix_method': 物理形态（'axial' | 'torsion' | 'bending_shear'）
%    - 'compare_analytical': 是否与解析解对比（默认 true）
%
%  保存图形（可选）：
%    visualization.saveFigure(h, 'filename')
%    visualization.saveFigure(h, 'filename', 'format', 'png')
%    visualization.saveFigure(h, 'filename', 'format', {'pdf', 'png'}, 'dpi', 300)
%    visualization.saveFigure(h, 'filename', 'format', 'pdf', 'path', 'output_dir')
%
%  saveFigure 参数说明：
%    - 'format': 文件格式（'pdf' | 'png' | 'eps' | 'svg' | 'fig'，默认'pdf'）
%    - 'dpi': PNG分辨率（默认300）
%    - 'path': 保存路径（默认当前目录）
% ============================================================

fprintf('\n======== 可视化示例 ========\n');

% 绘制弯曲-剪切模态（第1和第2阶）
% plotMode 会自动从缓存加载数据
h = visualization.plotMode('PP', beam, 100, [1 2], ...
    'matrix_method', 'bending_shear', ...
    'compare_analytical', true);

% 保存图形（可选）
visualization.saveFigure(h, 'PP_N100_mode_1_2_bending_shear', 'format', 'png');

fprintf('绘图完成！图形已保存为 PP_N100_mode_1_2_bending_shear.png\n');


%% ============================================================
%  附加说明：
%
%  1. 解析解的限制：
%     弯曲-剪切解析解在 λ≈α 附近存在数值奇异问题（α = AL²/I）。
%     - α ≥ 500（细长梁）：无影响
%     - α < 500（短粗梁）：可能缺失部分模态
%     详见：docs/analytical_solution_limitations.md
%
%  2. 缓存管理：
%     手动删除缓存：
%       rmdir('.cache', 's');  % 删除所有缓存
%
%  3. 更多示例：
%     - 不同边界条件对比
%     - 网格收敛性测试
%     - 截面参数影响分析
%     详见项目文档和 tests/ 目录
% ============================================================

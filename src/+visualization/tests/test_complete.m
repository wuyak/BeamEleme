%% test_complete_visualization - 可视化模块完整功能测试
%
% 目的：验证 visualization.plotMode() 的所有功能
% 测试方法：生成图像，手动视觉检查
%
% 检查项：
%   1. 所有物理类型都能正常绘图（axial, torsion, bending_shear）
%   2. FEM和解析解能够正确对比
%   3. 符号一致性（无异常翻转）
%   4. 多模态绘制
%   5. DOF选择功能
%   6. 边界条件支持

close all; clc;
fprintf('========================================\n');
fprintf('  可视化模块完整功能测试\n');
fprintf('========================================\n\n');

%% 配置
output_dir = fullfile(fileparts(mfilename('fullpath')), 'output');
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

NElem = 20;        % 单元数
n_modes = 10;      % 模态数

test_count = 0;
pass_count = 0;
fail_count = 0;

%% 准备梁模型
fprintf('准备梁模型...\n');
mat = parameters.MaterialLibrary.Steel();
sec_rect = parameters.SectionLibrary.Rectangular(0.2, 0.1);
beam = parameters.BeamModel(mat, sec_rect, 'L', 1.0);
fprintf('  ✓ 钢梁: L=1.0m, 0.2×0.1m矩形截面\n\n');

%% Test 1: Bending_Shear - PP边界 (多DOF)
fprintf('[Test 1] Bending_Shear - PP边界 (auto DOF)\n');
test_count = test_count + 1;
try
    cache = visualization.loadCache('PP', beam, NElem, n_modes);
    h1 = visualization.plotMode('PP', beam, NElem, [1 2 3], ...
        'matrix_method', 'bending_shear', 'dof', 'auto');
    saveas(h1, fullfile(output_dir, 'test1_bending_shear_PP_mode123_auto.png'));
    fprintf('  ✓ 图像已保存: test1_bending_shear_PP_mode123_auto.png\n');
    pass_count = pass_count + 1;
catch ME
    fprintf('  ✗ 失败: %s\n', ME.message);
    fail_count = fail_count + 1;
end

%% Test 2: Bending_Shear - PP边界 (只显示w)
fprintf('\n[Test 2] Bending_Shear - PP边界 (只显示w)\n');
test_count = test_count + 1;
try
    h2 = visualization.plotMode('PP', beam, NElem, [1 2], ...
        'matrix_method', 'bending_shear', 'dof', 'w');
    saveas(h2, fullfile(output_dir, 'test2_bending_shear_PP_mode12_w_only.png'));
    fprintf('  ✓ 图像已保存: test2_bending_shear_PP_mode12_w_only.png\n');
    pass_count = pass_count + 1;
catch ME
    fprintf('  ✗ 失败: %s\n', ME.message);
    fail_count = fail_count + 1;
end

%% Test 3: Bending_Shear - FF边界
fprintf('\n[Test 3] Bending_Shear - FF边界\n');
test_count = test_count + 1;
try
    cache = visualization.loadCache('FF', beam, NElem, n_modes);
    h3 = visualization.plotMode('FF', beam, NElem, [1 2 3], ...
        'matrix_method', 'bending_shear');
    saveas(h3, fullfile(output_dir, 'test3_bending_shear_FF_mode123.png'));
    fprintf('  ✓ 图像已保存: test3_bending_shear_FF_mode123.png\n');
    pass_count = pass_count + 1;
catch ME
    fprintf('  ✗ 失败: %s\n', ME.message);
    fail_count = fail_count + 1;
end

%% Test 4: Bending_Shear - CF边界 (悬臂梁)
fprintf('\n[Test 4] Bending_Shear - CF边界\n');
test_count = test_count + 1;
try
    cache = visualization.loadCache('CF', beam, NElem, n_modes);
    h4 = visualization.plotMode('CF', beam, NElem, [1 3 5], ...
        'matrix_method', 'bending_shear');
    saveas(h4, fullfile(output_dir, 'test4_bending_shear_CF_mode135.png'));
    fprintf('  ✓ 图像已保存: test4_bending_shear_CF_mode135.png\n');
    pass_count = pass_count + 1;
catch ME
    fprintf('  ✗ 失败: %s\n', ME.message);
    fail_count = fail_count + 1;
end

%% Test 5: Axial - PP边界
fprintf('\n[Test 5] Axial - PP边界\n');
test_count = test_count + 1;
try
    cache = visualization.loadCache('PP', beam, NElem, n_modes);
    h5 = visualization.plotMode('PP', beam, NElem, [1 3 5], 'matrix_method', 'axial');
    saveas(h5, fullfile(output_dir, 'test5_axial_PP_mode135.png'));
    fprintf('  ✓ 图像已保存: test5_axial_PP_mode135.png\n');
    pass_count = pass_count + 1;
catch ME
    fprintf('  ✗ 失败: %s\n', ME.message);
    fail_count = fail_count + 1;
end

%% Test 6: Axial - FF边界
fprintf('\n[Test 6] Axial - FF边界\n');
test_count = test_count + 1;
try
    cache = visualization.loadCache('FF', beam, NElem, n_modes);
    h6 = visualization.plotMode('FF', beam, NElem, [2 4], 'matrix_method', 'axial');
    saveas(h6, fullfile(output_dir, 'test6_axial_FF_mode24.png'));
    fprintf('  ✓ 图像已保存: test6_axial_FF_mode24.png\n');
    pass_count = pass_count + 1;
catch ME
    fprintf('  ✗ 失败: %s\n', ME.message);
    fail_count = fail_count + 1;
end

%% Test 7: Torsion - PP边界
fprintf('\n[Test 7] Torsion - PP边界\n');
test_count = test_count + 1;
try
    cache = visualization.loadCache('PP', beam, NElem, n_modes);
    h7 = visualization.plotMode('PP', beam, NElem, [1 2 3], 'matrix_method', 'torsion');
    saveas(h7, fullfile(output_dir, 'test7_torsion_PP_mode123.png'));
    fprintf('  ✓ 图像已保存: test7_torsion_PP_mode123.png\n');
    pass_count = pass_count + 1;
catch ME
    fprintf('  ✗ 失败: %s\n', ME.message);
    fail_count = fail_count + 1;
end

%% Test 8: Torsion - FF边界
fprintf('\n[Test 8] Torsion - FF边界\n');
test_count = test_count + 1;
try
    cache = visualization.loadCache('FF', beam, NElem, n_modes);
    h8 = visualization.plotMode('FF', beam, NElem, [2], 'matrix_method', 'torsion');
    saveas(h8, fullfile(output_dir, 'test8_torsion_FF_mode2.png'));
    fprintf('  ✓ 图像已保存: test8_torsion_FF_mode2.png\n');
    pass_count = pass_count + 1;
catch ME
    fprintf('  ✗ 失败: %s\n', ME.message);
    fail_count = fail_count + 1;
end

%% Test 9: 单个模态绘制
fprintf('\n[Test 9] 单个模态绘制 (Bending_Shear mode 1)\n');
test_count = test_count + 1;
try
    h9 = visualization.plotMode('PP', beam, NElem, 1, 'matrix_method', 'bending_shear');
    saveas(h9, fullfile(output_dir, 'test9_bending_shear_PP_mode1_single.png'));
    fprintf('  ✓ 图像已保存: test9_bending_shear_PP_mode1_single.png\n');
    pass_count = pass_count + 1;
catch ME
    fprintf('  ✗ 失败: %s\n', ME.message);
    fail_count = fail_count + 1;
end

%% Test 10: 不对比解析解 (仅FEM)
fprintf('\n[Test 10] 仅绘制FEM结果（不对比解析解）\n');
test_count = test_count + 1;
try
    h10 = visualization.plotMode('PP', beam, NElem, [1 2 3], ...
        'matrix_method', 'axial', 'compare_analytical', false);
    saveas(h10, fullfile(output_dir, 'test10_axial_PP_mode123_fem_only.png'));
    fprintf('  ✓ 图像已保存: test10_axial_PP_mode123_fem_only.png\n');
    pass_count = pass_count + 1;
catch ME
    fprintf('  ✗ 失败: %s\n', ME.message);
    fail_count = fail_count + 1;
end

%% 测试总结
fprintf('\n========================================\n');
fprintf('  测试总结\n');
fprintf('========================================\n');
fprintf('总测试数:   %d\n', test_count);
fprintf('通过:       %d (%.1f%%)\n', pass_count, pass_count/test_count*100);
fprintf('失败:       %d (%.1f%%)\n', fail_count, fail_count/test_count*100);
fprintf('\n生成的图像保存在:\n  %s\n', output_dir);
fprintf('\n请手动检查以下内容：\n');
fprintf('  1. ✓ 所有图像正常显示\n');
fprintf('  2. ✓ FEM和解析解曲线重合良好\n');
fprintf('  3. ✓ 符号一致（无异常翻转）\n');
fprintf('  4. ✓ 标题信息正确（模态号、频率、边界条件）\n');
fprintf('  5. ✓ 多DOF图像布局合理\n');
fprintf('========================================\n');

% 自动打开输出目录（Windows）
if ispc
    winopen(output_dir);
end

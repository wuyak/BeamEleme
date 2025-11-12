% demo_new_api - 重构后 visualization 模块快速演示
%
% 展示最常见的使用场景

clear; clc; close all;

% 添加src到路径
script_dir = fileparts(mfilename('fullpath'));         % src/+visualization/tests
vis_package_dir = fileparts(script_dir);               % src/+visualization
src_dir = fileparts(vis_package_dir);                  % src
addpath(src_dir);

fprintf('\n=== visualization API 快速演示 ===\n\n');

%% 1. 准备模型
fprintf('1. 准备梁模型...\n');
mat = parameters.MaterialLibrary.Steel();
sec = parameters.SectionLibrary.Rectangular(0.2, 0.1);
beam = parameters.BeamModel(mat, sec, 'L', 1.0);
fprintf('   ✓ 钢梁: L=1.0m, 0.2×0.1m矩形截面\n\n');

%% 2. 加载缓存
fprintf('2. 加载缓存（自动查找或生成）...\n');
cache = visualization.loadCache('PP', beam, 20, 3);
fprintf('   ✓ 缓存路径: %s\n', cache.cache_path);
fprintf('   ✓ 包含: axial, torsion, bending_shear\n\n');

%% 3. 查看数据
fprintf('3. 直接访问数据...\n');
fprintf('   弯曲剪切模态频率:\n');
for i = 1:3
    fem_freq = cache.bending_shear.fem.freq(i);
    ana_freq = cache.bending_shear.ana.freq_Hz(i);
    mac = cache.bending_shear.comp.mac(i);
    error = cache.bending_shear.comp.freq_error(i);
    fprintf('     模态%d: %.2f Hz (FEM) vs %.2f Hz (解析解)\n', i, fem_freq, ana_freq);
    fprintf('             误差=%.3f%%, MAC=%.4f\n', error, mac);
end
fprintf('\n');

%% 4. 绘制弯曲模态
fprintf('4. 绘制弯曲剪切第1阶模态...\n');
h1 = visualization.plotMode('PP', beam, 20, 1, 'matrix_method', 'bending_shear');
fprintf('   ✓ 自动显示2个子图: w(挠度) + theta(转角)\n');
fprintf('   ✓ 蓝色实线=FEM, 红色虚线=解析解\n');
fprintf('   ✓ 标题显示频率和MAC值\n\n');
pause(2);

%% 5. 绘制扭转模态
fprintf('5. 绘制扭转第1阶模态...\n');
h2 = visualization.plotMode('PP', beam, 20, 1, 'matrix_method', 'torsion');
fprintf('   ✓ 显示1个子图: phi(扭转角)\n\n');
pause(2);

%% 6. 自定义绘图
fprintf('6. 自定义绘图（只显示w DOF）...\n');
h3 = visualization.plotMode('PP', beam, 20, 2, 'matrix_method', 'bending_shear', ...
    'dof', 'w');
fprintf('   ✓ 只显示w位移（归一化显示）\n\n');
pause(2);

%% 7. 批量绘图示例
fprintf('7. 批量绘图示例（已跳过）\n');
fprintf('   提示: 使用循环绘制多个模态\n');
fprintf('   >> for i = 1:3\n');
fprintf('        h = visualization.plotMode(''bending_shear'', i, cache);\n');
fprintf('        saveFigure(h, sprintf(''mode_%%d.png'', i));\n');
fprintf('      end\n\n');

%% 总结
fprintf('=== 演示完成 ===\n');
fprintf('核心步骤:\n');
fprintf('  1. cache = loadCache(BC, beam, NElem, n_modes)\n');
fprintf('  2. h = plotMode(BC, beam, NElem, mode_idx, ''matrix_method'', type)\n');
fprintf('  3. 直接访问: cache.bending_shear.fem/ana/comp\n\n');

fprintf('当前打开的图形:\n');
fprintf('  - Figure %d: 弯曲剪切模态1 (w + theta)\n', h1.Number);
fprintf('  - Figure %d: 扭转模态1 (phi)\n', h2.Number);
fprintf('  - Figure %d: 弯曲剪切模态2 (仅w, 不归一化)\n', h3.Number);
fprintf('\n按任意键关闭所有图形...\n');

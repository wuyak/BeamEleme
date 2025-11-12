clc; clear;

% 添加项目路径
addpath(genpath('../src'));

%% 准备参数
mat = parameters.MaterialLibrary.Steel();
sec = parameters.SectionLibrary.Circular(0.6);
beam = parameters.BeamModel(mat, sec, 'L', 1.0);

%% 测试 1: 正常执行（完整报告）
fprintf('【测试 1】正常执行 - 完整报告\n');
fprintf('----------------------------------------\n');
try
    wf = workflow.solve('PP', beam, 500, 5);
    fprintf('✓ 测试通过\n\n');
catch ME
    fprintf('✗ 失败: %s\n\n', ME.message);
end

%% 测试 2: 输入验证 - 错误的 BC
fprintf('【测试 2】输入验证 - 错误的边界条件\n');
fprintf('----------------------------------------\n');
try
    wf = workflow.solve('XY', beam, 50, 5, false);
    fprintf('✗ 应该报错但没有报错\n\n');
catch ME
    fprintf('✓ 正确捕获错误: %s\n\n', ME.message);
end

%% 测试 3: 输入验证 - 负数 NElem
fprintf('【测试 3】输入验证 - 负数单元数\n');
fprintf('----------------------------------------\n');
try
    wf = workflow.solve('PP', beam, -5, 5, false);
    fprintf('✗ 应该报错但没有报错\n\n');
catch ME
    fprintf('✓ 正确捕获错误: %s\n\n', ME.message);
end

%% 测试 4: 输入验证 - 浮点数 n_modes
fprintf('【测试 4】输入验证 - 浮点数模态数\n');
fprintf('----------------------------------------\n');
try
    wf = workflow.solve('PP', beam, 50, 3.7, false);
    fprintf('✗ 应该报错但没有报错\n\n');
catch ME
    fprintf('✓ 正确捕获错误: %s\n\n', ME.message);
end

%% 测试 5: 输入验证 - 错误的 beam 类型
fprintf('【测试 5】输入验证 - 错误的 beam 类型\n');
fprintf('----------------------------------------\n');
try
    wf = workflow.solve('PP', struct('E', 2e11), 50, 5, false);
    fprintf('✗ 应该报错但没有报错\n\n');
catch ME
    fprintf('✓ 正确捕获错误: %s\n\n', ME.message);
end

%% 测试 6: 验证保存文件结构
fprintf('【测试 6】验证保存文件结构\n');
fprintf('----------------------------------------\n');
try
    wf = workflow.solve('CF', beam, 20, 3, false);

    % 检查保存路径
    assert(isfield(wf, 'save_path'), '缺少 save_path 字段');
    fprintf('  保存路径: %s\n', wf.save_path);

    % 验证目录存在
    assert(exist(wf.save_path, 'dir') == 7, '工作流目录不存在');

    % 验证子目录和文件
    subdirs = {'axial', 'torsion', 'bending_shear'};
    for i = 1:length(subdirs)
        subdir = fullfile(wf.save_path, subdirs{i});
        assert(exist(subdir, 'dir') == 7, sprintf('缺少 %s 子目录', subdirs{i}));

        % 验证文件
        fem_file = fullfile(subdir, 'fem_result.mat');
        ana_file = fullfile(subdir, 'ana_result.mat');
        assert(exist(fem_file, 'file') == 2, sprintf('缺少 %s/fem_result.mat', subdirs{i}));
        assert(exist(ana_file, 'file') == 2, sprintf('缺少 %s/ana_result.mat', subdirs{i}));
    end

    % 验证元数据文件
    meta_file = fullfile(wf.save_path, 'metadata.mat');
    assert(exist(meta_file, 'file') == 2, '缺少 metadata.mat');

    % 加载并验证元数据内容
    meta = load(meta_file);
    assert(isfield(meta, 'BC'), '缺少 metadata.BC');
    assert(isfield(meta, 'NElem'), '缺少 metadata.NElem');
    assert(isfield(meta, 'n_modes'), '缺少 metadata.n_modes');
    assert(strcmp(meta.BC, 'CF'), 'BC不匹配');
    assert(meta.NElem == 20, 'NElem不匹配');
    assert(meta.n_modes == 3, 'n_modes不匹配');

    % 测试加载单个物理类型结果
    fem_axial = load(fullfile(wf.save_path, 'axial', 'fem_result.mat'));
    assert(isfield(fem_axial, 'freq'), '缺少 FEM频率数据');
    assert(length(fem_axial.freq) == 3, 'FEM模态数不匹配');

    fprintf('✓ 保存文件结构正确\n');
    fprintf('  metadata.BC = %s\n', meta.BC);
    fprintf('  metadata.NElem = %d\n', meta.NElem);
    fprintf('  metadata.n_modes = %d\n', meta.n_modes);
    fprintf('  ✓ 所有子目录和文件完整\n\n');
catch ME
    fprintf('✗ 失败: %s\n\n', ME.message);
end

%% 总结
fprintf('========================================\n');
fprintf('测试完成！\n');
fprintf('========================================\n');

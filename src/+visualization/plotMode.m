function h = plotMode(BC, beam, NElem, mode_idx, varargin)
% plotMode - 从缓存读取并绘制指定模态振型
%
% 用法:
%   h = visualization.plotMode('CF', beam, 40, 6, 'matrix_method', 'bending_shear');
%   h = visualization.plotMode('PP', beam, 80, 3, 'matrix_method', 'axial');
%   h = visualization.plotMode('CC', beam, 40, [1 3 5], 'matrix_method', 'torsion');
%
% 输入:
%   BC         - 边界条件 ('PP', 'CF', 'CC', etc.)
%   beam       - parameters.BeamModel对象
%   NElem      - 单元数
%   mode_idx   - 模态索引（标量或向量，例如 6 或 [1 3 6]）
%
% Name-Value 参数:
%   'matrix_method' - 物理形态，默认 'bending_shear'
%                     可选: 'axial' | 'torsion' | 'bending_shear'
%   'cache_prefix'  - 缓存前缀，默认 ''
%   'dof'           - 指定DOF，默认 'auto'（自动选择）
%   'compare_analytical' - 是否与解析解对比，默认 false（实验性功能）
%
% 输出:
%   h - figure句柄
%
% 示例:
%   % 查看CF边界、弯曲剪切的第6个模态
%   mat = parameters.MaterialLibrary.Steel();
%   sec = parameters.SectionLibrary.Rectangular(0.2, 0.1);
%   beam = parameters.BeamModel(mat, sec, 'L', 1.0);
%   h = visualization.plotMode('CF', beam, 40, 6, 'matrix_method', 'bending_shear');
%
%   % 查看多个模态
%   h = visualization.plotMode('PP', beam, 40, [1 2 3], 'matrix_method', 'axial');

% 解析参数
p = inputParser;
addRequired(p, 'BC', @(x) ischar(x) || isstring(x));
addRequired(p, 'beam', @(x) isa(x, 'parameters.BeamModel'));
addRequired(p, 'NElem', @(x) isnumeric(x) && isscalar(x) && x > 0);
addRequired(p, 'mode_idx', @isnumeric);
addParameter(p, 'matrix_method', 'bending_shear', @(x) ischar(x) || isstring(x));
addParameter(p, 'cache_prefix', '', @(x) ischar(x) || isstring(x));
addParameter(p, 'dof', 'auto', @(x) ischar(x) || isstring(x));
addParameter(p, 'compare_analytical', true, @islogical);
parse(p, BC, beam, NElem, mode_idx, varargin{:});
% === 配置常量 ===
MAC_THRESHOLD = 0.8;           % 伪解判定阈值
SAMPLE_POINTS = 100;           % 解析解采样点数
FIG_WIDTH_PER_MODE = 300;      % 每个模态的图宽(像素)
FIG_HEIGHT_SINGLE = 400;       % 单DOF图高(像素)


BC = upper(char(p.Results.BC));
matrix_method = lower(char(p.Results.matrix_method));
cache_prefix = char(p.Results.cache_prefix);
dof = p.Results.dof;
compare_analytical = p.Results.compare_analytical;

% 验证 mode_idx
max_mode = max(mode_idx);
if any(mode_idx < 1)
    error('mode_idx 必须 >= 1');
end

% Step 1: 智能查找缓存（需要 M >= max(mode_idx)）
cache_path = findBestCache(BC, beam, NElem, max_mode, cache_prefix);

if isempty(cache_path)
    error('未找到满足 M>=%d 的缓存。请先运行:\n  workflow.solve(''%s'', beam, %d, %d)', ...
        max_mode, BC, NElem, max_mode);
end

% 显示缓存信息
[~, cache_name] = fileparts(cache_path);
fprintf('✓ 使用缓存: %s\n', cache_name);

% Step 2: 加载数据
fem_file = fullfile(cache_path, matrix_method, 'fem_result.mat');
if ~exist(fem_file, 'file')
    error('缓存中缺少 %s 数据。\n可用的物理形态: %s', ...
        matrix_method, strjoin(getFolders(cache_path), ', '));
end

fem_result = load(fem_file);

% 验证模态数量
n_modes_available = size(fem_result.modes, 2);
if max_mode > n_modes_available
    error('缓存只有 %d 个模态，无法提取第 %d 个', n_modes_available, max_mode);
end

% Step 3: 加载解析解（如果需要对比）
ana_result = [];
matching_info = [];
if compare_analytical
    ana_file = fullfile(cache_path, matrix_method, 'ana_result.mat');
    if exist(ana_file, 'file')
        ana_result = load(ana_file);
        fprintf('✓ 加载解析解数据\n');

        % 创建简化的matching_info（跳过MAC计算）
        for i = 1:length(mode_idx)
            matching_info(i).fem_idx = mode_idx(i);
            matching_info(i).ana_idx = mode_idx(i);  % 假设1:1对应
            matching_info(i).mac = 1.0;
            matching_info(i).freq_error = 0;
            matching_info(i).is_spurious = false;
            matching_info(i).is_missed = false;
        end
    else
        fprintf('⚠ 未找到解析解数据，仅绘制FEM结果\n');
        compare_analytical = false;
    end
end

% Step 4: 绘制模态振型
if compare_analytical && ~isempty(ana_result)
    % 绘制对比图
    h = plotModeComparison(fem_result, ana_result, mode_idx, dof, matching_info);
else
    % 仅绘制FEM
    h = plotFemOnly(fem_result, mode_idx, dof);
end

% 更新标题
mode_str = sprintf('%d', mode_idx(1));
if length(mode_idx) > 1
    mode_str = sprintf('[%s]', num2str(mode_idx));
end

% 生成梁特征描述
beam_desc = generateBeamDescription(beam, NElem);

title_str = sprintf('%s - %s 模态振型\n模态: %s | %s', ...
    BC, strrep(matrix_method, '_', '-'), mode_str, beam_desc);

% 添加匹配信息到标题
if compare_analytical && ~isempty(matching_info)
    for i = 1:length(mode_idx)
        if matching_info(i).is_spurious
            title_str = sprintf('%s\n模态%d: 数值伪解 (MAC<0.8)', mode_idx(i));
        elseif matching_info(i).mac < 0.95
            title_str = sprintf('%s\n模态%d: 匹配度%.2f%%', mode_idx(i), matching_info(i).mac*100);
        end
    end
end

sgtitle(gcf, title_str, 'FontSize', 12, 'FontWeight', 'bold');

fprintf('✓ 模态振型绘制完成。\n');

end


function folders = getFolders(cache_path)
% getFolders - 获取缓存中可用的物理形态
d = dir(cache_path);
folders = {d([d.isdir] & ~strncmp({d.name}, '.', 1)).name};
end


function matching_info = analyzeModeMatching(fem_result, ana_result, mode_idx)
% analyzeModeMatching - 分析FEM和解析解的模态匹配情况
%
% 输出:
%   matching_info(i).fem_idx        - FEM模态索引
%   matching_info(i).ana_idx        - 匹配的解析解索引（0表示无匹配）
%   matching_info(i).mac            - MAC值
%   matching_info(i).freq_error     - 频率误差 (%)
%   matching_info(i).is_spurious    - 是否为数值伪解
%   matching_info(i).is_missed      - 对应的解析解模态是否被FEM漏掉


n_fem = size(fem_result.modes, 2);
n_ana = size(ana_result.modes, 2);

% 计算完整的MAC矩阵
MAC = zeros(n_fem, n_ana);
for i = 1:n_fem
    for j = 1:n_ana
        % FEM: 直接提取列
        fem_mode = fem_result.modes(:,i);

        % 解析解: 可能是 cell{1,j}{dof_idx,1}
        if iscell(ana_result.modes)
            ana_mode = ana_result.modes{1,j}{1,1};  % 假设第一个DOF
        else
            ana_mode = ana_result.modes(:,j);
        end

        MAC(i,j) = computeMAC(fem_mode, ana_mode);
    end
end

% 分析每个请求的模态
matching_info = struct('fem_idx', {}, 'ana_idx', {}, 'mac', {}, ...
    'freq_error', {}, 'is_spurious', {}, 'is_missed', {});

for i = 1:length(mode_idx)
    idx = mode_idx(i);

    if idx > n_fem
        continue;
    end

    % 找到最佳匹配的解析解模态
    [max_mac, best_ana] = max(MAC(idx, :));

    % 计算频率误差
    freq_fem = fem_result.freq(idx);
    freq_ana = ana_result.freq(best_ana);
    freq_error = abs(freq_fem - freq_ana) / freq_ana * 100;

    % 判断是否为伪解
    is_spurious = (max_mac < MAC_THRESHOLD);

    matching_info(i).fem_idx = idx;
    matching_info(i).ana_idx = best_ana;
    matching_info(i).mac = max_mac;
    matching_info(i).freq_error = freq_error;
    matching_info(i).is_spurious = is_spurious;
    matching_info(i).is_missed = false;

    % 输出诊断信息
    if is_spurious
        fprintf('  ⚠ 模态%d: 数值伪解 (MAC=%.3f < %.2f)\n', idx, max_mac, MAC_THRESHOLD);
    elseif max_mac < 0.95
        fprintf('  ⚠ 模态%d: 匹配度偏低 (MAC=%.3f, 频率误差=%.2f%%)\n', ...
            idx, max_mac, freq_error);
    else
        fprintf('  ✓ 模态%d: 匹配良好 (MAC=%.4f, 频率误差=%.2f%%)\n', ...
            idx, max_mac, freq_error);
    end
end

% 检查是否有解析解模态被漏掉（miss）
matched_ana = [matching_info.ana_idx];
for j = 1:min(n_ana, max(mode_idx))
    if ~ismember(j, matched_ana)
        fprintf('  ⚠ 解析解模态%d可能被FEM漏掉 (miss)\n', j);
    end
end
end


function mac = computeMAC(mode1, mode2)
% computeMAC - 计算模态保证准则
%
% MAC = |φ1' * φ2|² / [(φ1' * φ1) * (φ2' * φ2)]

% 统一采样点
n_sample = SAMPLE_POINTS;
x = linspace(0, 1, n_sample)';

% 采样 mode1
if isa(mode1, 'chebfun')
    v1 = feval(mode1, x);
    if size(v1, 2) > 1  % 多列（多DOF），取第一列
        v1 = v1(:, 1);
    end
elseif isnumeric(mode1)
    v1 = mode1(:);
    if length(v1) ~= n_sample
        x1 = linspace(0, 1, length(v1))';
        v1 = interp1(x1, v1, x, 'linear');
    end
else
    error('Unsupported mode1 type: %s', class(mode1));
end

% 采样 mode2
if isa(mode2, 'chebfun')
    v2 = feval(mode2, x);
    if size(v2, 2) > 1  % 多列（多DOF），取第一列
        v2 = v2(:, 1);
    end
elseif isnumeric(mode2)
    v2 = mode2(:);
    if length(v2) ~= n_sample
        x2 = linspace(0, 1, length(v2))';
        v2 = interp1(x2, v2, x, 'linear');
    end
else
    error('Unsupported mode2 type: %s', class(mode2));
end

% 确保是列向量
v1 = v1(:);
v2 = v2(:);

% 计算MAC
numerator = abs(v1' * v2)^2;
denominator = (v1' * v1) * (v2' * v2);

if denominator == 0
    mac = 0;
else
    mac = numerator / denominator;
end
end


function h = plotModeComparison(fem_result, ana_result, mode_idx, dof, matching_info)
FIG_WIDTH_PER_MODE = 300;
FIG_HEIGHT_SINGLE = 400;
% plotModeComparison - 绘制FEM和解析解对比图

n_modes = length(mode_idx);

% 确定DOF列表
if strcmp(dof, 'auto')
    % 根据物理类型确定默认DOF列表
    if strcmp(fem_result.matrix_method, 'bending_shear')
        dof_list = {'w', 'theta'};
    elseif strcmp(fem_result.matrix_method, 'axial')
        dof_list = {'u'};
    else  % torsion
        dof_list = {'phi'};
    end
else
    dof_list = {dof};
end

n_dofs = length(dof_list);

% 创建figure（bending_shear有2个DOF，需要2行）
if n_dofs == 1
    h = figure('Position', [100 100 FIG_WIDTH_PER_MODE*n_modes FIG_HEIGHT_SINGLE]);
else
    h = figure('Position', [100 100 FIG_WIDTH_PER_MODE*n_modes FIG_HEIGHT_SINGLE*n_dofs]);
end

for i = 1:n_modes
    idx = mode_idx(i);

    for j = 1:n_dofs
        dof_name = dof_list{j};

        % 子图索引
        if n_dofs == 1
            subplot(1, n_modes, i);
        else
            subplot(n_dofs, n_modes, (j-1)*n_modes + i);
        end

        hold on; grid on;

        % 提取FEM数据
        [x_fem, y_fem] = extractDOF(fem_result, dof_name, idx);

        % 提取解析解数据
        [x_ana, y_ana] = extractAnalyticalMode(ana_result, dof_name, idx, matching_info(i));

        % 总是归一化以比较振型形状（不论normalize参数）
        % 振型比较关注形状，不关注幅值
        norm_fem = max(abs(y_fem));
        if norm_fem > 0
            y_fem = y_fem / norm_fem;
        end

        if ~isempty(y_ana)
            norm_ana = max(abs(y_ana));
            if norm_ana > 0
                y_ana = y_ana / norm_ana;
            end
        end

        % 符号调整：确保 FEM 和解析解符号一致
        % 简单粗暴：插值解析解到FEM节点，计算相关性，小于0就翻转
        if ~isempty(y_ana) && length(x_ana) > 1 && length(x_fem) > 1
            y_ana_sampled = interp1(x_ana, y_ana, x_fem, 'linear', 'extrap');
            correlation = y_fem' * y_ana_sampled;
            if correlation < 0
                y_ana = -y_ana;  % 翻转整个解析解
            end
        end

        % 检查解析解有效性
        has_ana = ~isempty(x_ana) && ~isempty(y_ana) && isnumeric(y_ana) && ~all(isnan(y_ana));

        % 绘制并显式捕获句柄
        h_fem = plot(x_fem, y_fem, 'b-', 'LineWidth', 2);

        if has_ana
            % 使用粗虚线，确保可见
            h_ana = plot(x_ana, y_ana, 'r--', 'LineWidth', 3);
            legend([h_fem h_ana], {'FEM', '解析解'}, 'Location', 'best');
        else
            legend(h_fem, {'FEM'}, 'Location', 'best');
        end

        % 标题和标签
        freq_fem = fem_result.freq(idx);
        if j == 1
            % 第一个DOF显示模态号和频率
            title_text = sprintf('模态 %d - %.1f Hz', idx, freq_fem);
        else
            % 其他DOF只显示DOF名称
            title_text = sprintf('%s', getDOFLabel(dof_name));
        end

        % 添加匹配信息（只在第一个DOF显示）
        if j == 1 && ~isempty(matching_info) && i <= length(matching_info)
            if matching_info(i).is_spurious
                title_text = sprintf('%s\n数值伪解 (MAC=%.2f)', title_text, matching_info(i).mac);
            elseif matching_info(i).mac < 0.95
                title_text = sprintf('%s\nMAC=%.3f', title_text, matching_info(i).mac);
            end
        end

        title(title_text);
        xlabel('位置 [m]');
        ylabel(getDOFLabel(dof_name));

        hold off;
    end
end
end


function [x, y] = extractAnalyticalMode(ana_result, dof_name, fem_idx, match_info)
% extractAnalyticalMode - 提取对应的解析解模态（简化版本）
%
% 重构后的简洁实现：3层嵌套，逻辑清晰

% 第1层：数据有效性检查
if isempty(ana_result) || match_info.is_spurious
    x = [];
    y = [];
    return;
end

ana_idx = match_info.ana_idx;

% 获取坐标
if isfield(ana_result, 'params') && isfield(ana_result.params, 'L')
    L = ana_result.params.L;
else
    L = 1.0;
end

try
    % 第2层：根据数据结构类型分发
    if ~isfield(ana_result, 'modes')
        % Case 1: 使用getModes函数（Axial/Torsion）
        if isfield(ana_result, 'getModes')
            [x, y] = handleGetModesCase(ana_result, dof_name, ana_idx, L);
        else
            x = [];
            y = [];
        end
        
    elseif iscell(ana_result.modes)
        % 检测是否为bending_shear的multi-region结构
        is_multi_region = (length(ana_result.modes) == 3) && ...
                          isfield(ana_result, 'lambda') && ...
                          isfield(ana_result, 'sys');
        
        if is_multi_region
            % Case 2: Multi-region结构（Bending_Shear）
            [x, y] = handleMultiRegion(ana_result, dof_name, ana_idx, L);
        else
            % Case 3: 单region结构
            [x, y] = handleSingleRegion(ana_result, dof_name, ana_idx, L);
        end
        
    else
        % Case 4: 非cell结构（直接数组）
        x = linspace(0, L, 100)';
        % 计算DOF索引
        if strcmp(ana_result.matrix_method, 'bending_shear')
            if strcmp(lower(dof_name), 'w')
                indices = 1:2:size(ana_result.modes, 1);  % 奇数行
            else
                indices = 2:2:size(ana_result.modes, 1);  % 偶数行
            end
        else
            indices = 1:size(ana_result.modes, 1);
        end
        y = ana_result.modes(indices(1), ana_idx);
        if length(y) ~= length(x)
            x_old = linspace(0, L, length(y))';
            y = interp1(x_old, y(:), x, 'linear');
        end
        y = y(:);
    end
    
catch ME
    % 第3层：错误处理
    fprintf('  ⚠ 解析解提取失败 (%s, mode %d): %s\n', dof_name, ana_idx, ME.message);
    x = [];
    y = [];
end
end
function desc = generateBeamDescription(beam, NElem)
% generateBeamDescription - 生成梁的特征描述
%
% 输出示例:
%   "Steel | 0.2×0.1m | L=1.0m | N=40"
%   "Aluminum | Circular D=0.05m | L=2.5m | N=80"

% 材料名称（简化）
material_name = getMaterialName(beam);

% 截面描述
section_desc = getSectionDescription(beam);

% 组合描述
desc = sprintf('%s | %s | L=%.1fm | N=%d', ...
    material_name, section_desc, beam.L, NElem);
end


function name = getMaterialName(beam)
% getMaterialName - 从BeamModel获取材料名称
%
% 根据常见材料属性判断

if beam.E == 2e11
    name = 'Steel';
elseif beam.E == 7e10
    name = 'Aluminum';
elseif beam.E == 1e10
    name = 'Concrete';
else
    % 默认显示弹性模量
    name = sprintf('E=%.1eN/m²', beam.E);
end
end


function desc = getSectionDescription(beam)
% getSectionDescription - 从BeamModel获取截面描述

% 截面类型首字母大写
section_name = beam.section_type;
if ~isempty(section_name)
    section_name(1) = upper(section_name(1));
end

if strcmp(beam.section_type, 'rectangular')
    % 矩形截面：从面积和惯性矩反推尺寸
    % A = b*h, I = b*h³/12
    if isfield(beam.dims, 'b') && isfield(beam.dims, 'h')
        desc = sprintf('%s %.2g×%.2gm', section_name, beam.dims.b, beam.dims.h);
    else
        % 从 A 和 I 估算
        h = (12 * beam.I / beam.A)^(1/2);
        b = beam.A / h;
        desc = sprintf('%s %.2g×%.2gm', section_name, b, h);
    end

elseif strcmp(beam.section_type, 'circular')
    % 圆形截面：A = π*r², I = π*r⁴/4
    if isfield(beam.dims, 'r')
        desc = sprintf('%s D=%.2gm', section_name, 2*beam.dims.r);
    else
        r = (beam.I * 4 / pi)^(1/4);
        desc = sprintf('%s D=%.2gm', section_name, 2*r);
    end

else
    % 其他截面类型，显示面积
    desc = sprintf('%s A=%.2gm²', section_name, beam.A);
end
end

function h = plotFemOnly(fem_result, mode_idx, dof)
FIG_WIDTH_PER_MODE = 300;
FIG_HEIGHT_SINGLE = 400;
% plotFemOnly - 仅绘制FEM模态（无解析解对比）
%
% 这是 plotModeComparison 的简化版本，用于没有解析解的情况

n_modes = length(mode_idx);

% 确定DOF列表
if strcmp(dof, 'auto')
    % 根据物理类型确定默认DOF列表
    if strcmp(fem_result.matrix_method, 'bending_shear')
        dof_list = {'w', 'theta'};
    elseif strcmp(fem_result.matrix_method, 'axial')
        dof_list = {'u'};
    else  % torsion
        dof_list = {'phi'};
    end
else
    dof_list = {dof};
end

n_dofs = length(dof_list);

% 创建figure
if n_dofs == 1
    h = figure('Position', [100 100 FIG_WIDTH_PER_MODE*n_modes FIG_HEIGHT_SINGLE]);
else
    h = figure('Position', [100 100 FIG_WIDTH_PER_MODE*n_modes FIG_HEIGHT_SINGLE*n_dofs]);
end

for i = 1:n_modes
    idx = mode_idx(i);
    
    for j = 1:n_dofs
        dof_name = dof_list{j};
        
        % 子图索引
        if n_dofs == 1
            subplot(1, n_modes, i);
        else
            subplot(n_dofs, n_modes, (j-1)*n_modes + i);
        end
        
        hold on; grid on;
        
        % 提取FEM数据
        [x_fem, y_fem] = extractDOF(fem_result, dof_name, idx);

        % 归一化（模态振型总是归一化以便比较）
        norm_fem = max(abs(y_fem));
        if norm_fem > 0
            y_fem = y_fem / norm_fem;
        end

        % 绘制FEM
        plot(x_fem, y_fem, 'b-', 'LineWidth', 2);
        legend({'FEM'}, 'Location', 'best');
        
        % 标题和标签
        freq_fem = fem_result.freq(idx);
        if j == 1
            title_text = sprintf('模态 %d - %.1f Hz', idx, freq_fem);
        else
            title_text = sprintf('%s', getDOFLabel(dof_name));
        end
        
        title(title_text);
        xlabel('位置 [m]');
        ylabel(getDOFLabel(dof_name));
        
        hold off;
    end
end
end


function cache_path = findBestCache(BC, beam, NElem, n_modes, cache_prefix)
% findBestCache - 智能查找满足 M>=n_modes 的最佳缓存
%
% 查找策略：
%   1. 精确匹配：M==n_modes
%   2. 最小满足：M>=n_modes 中最小的M
%   3. 未找到：返回空

% 获取项目根目录
script_path = mfilename('fullpath');
viz_dir = fileparts(script_path);       % +visualization/
src_dir = fileparts(viz_dir);           % src/
project_root = fileparts(src_dir);      % project root/

% 构建搜索路径
param_hash = beam.toHash();
search_pattern = sprintf('%s_%s_N%d_M*', param_hash(1:8), BC, NElem);

if isempty(cache_prefix)
    search_dir = fullfile(project_root, '.cache');
else
    search_dir = fullfile(project_root, '.cache', cache_prefix);
end

% 查找所有匹配的缓存目录
if ~exist(search_dir, 'dir')
    cache_path = '';
    return;
end

all_dirs = dir(fullfile(search_dir, search_pattern));
all_dirs = all_dirs([all_dirs.isdir]);  % 只保留目录

if isempty(all_dirs)
    cache_path = '';
    return;
end

% 解析模态数并筛选
candidates = struct('path', {}, 'M', {});
for i = 1:length(all_dirs)
    % 从目录名提取M值：{hash}_{BC}_N{NElem}_M{M}
    tokens = regexp(all_dirs(i).name, '_M(\d+)$', 'tokens');
    if ~isempty(tokens)
        M = str2double(tokens{1}{1});
        if M >= n_modes
            candidates(end+1).path = fullfile(search_dir, all_dirs(i).name);
            candidates(end).M = M;
        end
    end
end

% 未找到满足条件的缓存
if isempty(candidates)
    cache_path = '';
    return;
end

% 排序：选择M最小的（最节省内存，最快加载）
[~, idx] = min([candidates.M]);
cache_path = candidates(idx).path;

end


function [x, y] = handleGetModesCase(ana_result, dof_name, ana_idx, L)
% handleGetModesCase - 处理使用getModes函数的情况（Axial/Torsion）

x = linspace(0, L, 100)';
modes_sampled = ana_result.getModes(x / L);  % 归一化坐标

if size(modes_sampled, 2) >= ana_idx
    y = modes_sampled(:, ana_idx);
    y = y(:);
else
    x = [];
    y = [];
end
end


function [x, y] = handleMultiRegion(ana_result, dof_name, ana_idx, L)
% handleMultiRegion - 处理bending_shear的3区域结构

x = linspace(0, L, 100)';

% 使用 packLambdaForModes 获取每个region的模态分布
lambda_cell = analytical.bending_shear.packLambdaForModes(...
    ana_result.lambda, ana_result.sys.par.alpha);
n_modes_per_region = cellfun(@length, lambda_cell);

% 根据 ana_idx 定位所属的region
cumsum_modes = [0, cumsum(n_modes_per_region(:)')];
region_idx = find(ana_idx > cumsum_modes, 1, 'last');
local_idx = ana_idx - cumsum_modes(region_idx);

% 确定DOF索引
if strcmp(dof_name, 'w') || strcmp(dof_name, 'u')
    dof_idx = 1;
elseif strcmp(dof_name, 'theta') || strcmp(dof_name, 'phi')
    dof_idx = 2;
else
    dof_idx = 1;
end

% 从对应region提取chebfun
region_modes = ana_result.modes{region_idx};
if ~iscell(region_modes) || length(region_modes) < dof_idx
    error('Region %d 缺少 DOF %d 的数据', region_idx, dof_idx);
end

cheb_obj = region_modes{dof_idx};

% 采样
if isa(cheb_obj, 'chebfun')
    vals = feval(cheb_obj, x / L);
    if size(vals, 2) >= local_idx
        y = vals(:, local_idx);
    else
        error('Region %d 的模态数不足（需要 %d，只有 %d）', ...
            region_idx, local_idx, size(vals, 2));
    end
else
    error('Region %d 的 DOF %d 不是 chebfun', region_idx, dof_idx);
end

y = y(:);
end


function [x, y] = handleSingleRegion(ana_result, dof_name, ana_idx, L)
% handleSingleRegion - 处理单region结构（普通cell结构）

x = linspace(0, L, 100)';

mode_cell = ana_result.modes{1, 1};

% 确定DOF索引
if strcmp(dof_name, 'w')
    dof_idx = 1;
elseif strcmp(dof_name, 'theta')
    dof_idx = 2;
elseif strcmp(dof_name, 'u')
    dof_idx = 1;
elseif strcmp(dof_name, 'phi')
    dof_idx = 1;
else
    dof_idx = 1;
end

% 第二层：处理cell或直接chebfun
if iscell(mode_cell)
    if length(mode_cell) < dof_idx
        error('dof_idx %d 超出范围', dof_idx);
    end
    cheb_obj = mode_cell{1, dof_idx};
else
    % mode_cell直接是chebfun
    cheb_obj = mode_cell;
end

% 采样
if isa(cheb_obj, 'chebfun')
    vals = feval(cheb_obj, x / L);
    if size(vals, 2) >= ana_idx
        y = vals(:, ana_idx);
    else
        y = vals(:, 1);  % 降级：使用第一列
    end
else
    y = cheb_obj(:);
end

y = y(:);
end


function [x, value] = extractDOF(result, dof_name, mode_idx)
% extractDOF - 提取特定DOF的模态数据（局部辅助函数）
%
% 替代 extractDOF，逻辑直接实现，无过度抽象

% 计算节点数和坐标
NNodes = result.NElem + 1;
x = result.params.L * (0:NNodes-1)' / (NNodes - 1);

% 提取DOF值
if strcmp(result.matrix_method, 'bending_shear')
    % 弯曲剪切：2-DOF交错存储 [w1, θ1, w2, θ2, ...]
    if strcmp(lower(dof_name), 'w')
        value = result.modes(1:2:end, mode_idx);  % 奇数行：w
    else  % theta
        value = result.modes(2:2:end, mode_idx);  % 偶数行：θ
    end
else
    % 轴向/扭转：1-DOF
    value = result.modes(:, mode_idx);
end
end


function label = getDOFLabel(dof_name)
% getDOFLabel - 生成DOF的ylabel标签（局部辅助函数）
%
% 替代 getDOFLabel，简单映射，无过度抽象

switch lower(dof_name)
    case 'u'
        label = 'u [m]';
    case 'w'
        label = 'w [m]';
    case {'theta', 'θ'}
        label = '\theta [rad]';
    case {'phi', 'φ'}
        label = '\phi [rad]';
    otherwise
        label = sprintf('%s', dof_name);
end
end

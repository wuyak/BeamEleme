function cache_data = loadCache(BC, beam, NElem, n_modes, varargin)
% loadCache - 统一缓存管理入口
%
% 功能：
%   1. 智能查找 M>=n_modes 的最小缓存
%   2. 若不存在 → 调用 workflow.solve() 生成
%   3. 加载三个物理类型的 fem/ana/comp 结果
%   4. 返回原始数据（不截取，绘图时按需索引）
%
% 输入：
%   BC           - 边界条件 ('PP', 'CC', 'CF', etc.)
%   beam         - BeamModel对象
%   NElem        - 单元数
%   n_modes      - 需要的模态数
%
% Name-Value参数：
%   'cache_prefix'    - 缓存子目录，默认 ''
%   'force_recompute' - 强制重新计算，默认 false
%   'matrix_methods'  - 要加载的物理类型，默认 'all'
%                       'all' | {'axial', 'torsion', 'bending_shear'}
%
% 输出：
%   cache_data.axial.fem          - 轴向FEM结果（完整M个模态）
%   cache_data.axial.ana          - 轴向解析解
%   cache_data.axial.comp         - 轴向对比结果
%   cache_data.torsion.fem/ana/comp
%   cache_data.bending_shear.fem/ana/comp
%   cache_data.metadata           - 元数据
%   cache_data.cache_path         - 缓存路径
%
% 示例：
%   mat = parameters.MaterialLibrary.Steel();
%   sec = parameters.SectionLibrary.Rectangular(0.2, 0.1);
%   beam = parameters.BeamModel(mat, sec, 'L', 1.0);
%
%   % 加载缓存（若不存在自动生成）
%   cache = visualization.loadCache('PP', beam, 40, 10);
%
%   % 访问数据
%   fem = cache.bending_shear.fem;
%   freq_3 = fem.freq(3);           % 第3阶频率
%   mode_3 = fem.modes(:, 3);       % 第3阶模态
%   mac_3 = cache.bending_shear.comp.mac(3);  % MAC值

% 解析参数
p = inputParser;
addRequired(p, 'BC', @(x) ischar(x) || isstring(x));
addRequired(p, 'beam', @(x) isa(x, 'parameters.BeamModel'));
addRequired(p, 'NElem', @(x) isnumeric(x) && isscalar(x) && x > 0);
addRequired(p, 'n_modes', @(x) isnumeric(x) && isscalar(x) && x > 0);
addParameter(p, 'cache_prefix', '', @(x) ischar(x) || isstring(x));
addParameter(p, 'force_recompute', false, @islogical);
addParameter(p, 'matrix_methods', 'all');
parse(p, BC, beam, NElem, n_modes, varargin{:});

BC = upper(char(p.Results.BC));
cache_prefix = char(p.Results.cache_prefix);
force_recompute = p.Results.force_recompute;
matrix_methods = p.Results.matrix_methods;

% 标准化 matrix_methods
if ischar(matrix_methods) || isstring(matrix_methods)
    if strcmpi(matrix_methods, 'all')
        matrix_methods = {'axial', 'torsion', 'bending_shear'};
    else
        matrix_methods = {char(matrix_methods)};
    end
end

%% 查找或生成缓存

if force_recompute
    cache_path = '';
else
    % 智能查找 M>=n_modes 的最小缓存
    cache_path = findBestCache(BC, beam, NElem, n_modes, cache_prefix);
end

if isempty(cache_path)
    % 未找到缓存，调用 workflow.solve 生成
    fprintf('未找到缓存，生成新缓存: BC=%s, NElem=%d, n_modes=%d\n', BC, NElem, n_modes);
    wf_result = workflow.solve(BC, beam, NElem, n_modes, false, cache_prefix);
    cache_path = wf_result.save_path;
else
    fprintf('✓ 找到缓存: %s\n', cache_path);
end

%% 加载缓存数据

% 加载元数据
meta_file = fullfile(cache_path, 'metadata.mat');
if exist(meta_file, 'file') ~= 2
    error('visualization:loadCache:MissingMetadata', ...
        '缓存缺少 metadata.mat: %s', cache_path);
end
metadata = load(meta_file);

% 初始化输出
cache_data.cache_path = cache_path;
cache_data.metadata = metadata;

% 加载每个物理形态
for i = 1:numel(matrix_methods)
    method = matrix_methods{i};
    method_dir = fullfile(cache_path, method);

    if exist(method_dir, 'dir') ~= 7
        warning('visualization:loadCache:MissingMethod', ...
            '缓存缺少 %s 子目录，跳过', method);
        cache_data.(method) = [];
        continue;
    end

    fem_file = fullfile(method_dir, 'fem_result.mat');
    ana_file = fullfile(method_dir, 'ana_result.mat');
    comp_file = fullfile(method_dir, 'comp_result.mat');

    % 加载 FEM 结果
    if exist(fem_file, 'file') ~= 2
        warning('visualization:loadCache:MissingFEM', ...
            '缺少 FEM 结果: %s', fem_file);
        cache_data.(method) = [];
        continue;
    end
    fem = load(fem_file);

    % 加载解析解
    if exist(ana_file, 'file') ~= 2
        warning('visualization:loadCache:MissingAnalytical', ...
            '缺少解析解: %s', ana_file);
        ana = struct();
    else
        ana = load(ana_file);
    end

    % 加载对比结果（优先从缓存加载）
    if exist(comp_file, 'file') == 2
        % 直接加载已保存的 comp
        comp = load(comp_file);
    else
        % 旧缓存没有 comp_result.mat，动态计算
        if ~isempty(ana)
            % 计算对比时使用缓存中的完整模态数
            cached_n_modes = metadata.n_modes;
            comp = comparison.compare(fem, ana, 'n_modes', cached_n_modes, 'verbose', false);
        else
            comp = struct();
        end
    end

    % 打包（保留完整数据，不截取）
    cache_data.(method).fem = fem;
    cache_data.(method).ana = ana;
    cache_data.(method).comp = comp;
end

end


%% 辅助函数：查找最佳缓存

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

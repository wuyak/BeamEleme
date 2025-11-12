function wf_result = solve(BC, beam, NElem, n_modes, show_report, cache_prefix)
% solve - 完整求解流程：FEM → 解析解 → 对比 → 保存（所有物理形态）
%
% 用法:
%   wf_result = workflow.solve(BC, beam, NElem, n_modes)
%   wf_result = workflow.solve(BC, beam, NElem, n_modes, show_report)
%   wf_result = workflow.solve(BC, beam, NElem, n_modes, show_report, cache_prefix)
%
% 输入:
%   BC            - 边界条件字符串 ('PP', 'CC', 'CF', 'FF', etc.)
%   beam          - BeamModel对象（只包含物理属性）
%   NElem         - 单元数
%   n_modes       - 模态数
%   show_report   - (可选) 是否输出报告，默认 true
%   cache_prefix  - (可选) 缓存文件夹名称，默认为空
%                   留空: .cache/{hash}_{BC}_N{NElem}_M{n_modes}/
%                   指定: .cache/{cache_prefix}/{hash}_{BC}_N{NElem}_M{n_modes}/
%
% 输出:
%   wf_result.axial         - 轴向振动结果 (fem, ana, comp)
%   wf_result.torsion       - 扭转振动结果 (fem, ana, comp)
%   wf_result.bending_shear - 弯曲剪切振动结果 (fem, ana, comp)
%   wf_result.save_path     - 保存路径
%
% 示例:
%   mat = parameters.MaterialLibrary.Steel();
%   sec = parameters.SectionLibrary.Rectangular(0.2, 0.1);
%   beam = parameters.BeamModel(mat, sec, 'L', 1.0);
%
%   % 默认输出报告，缓存在 .cache/ 根目录
%   wf_result = workflow.solve('PP', beam, 100, 10);
%
%   % 静默模式（批量执行）
%   wf_result = workflow.solve('PP', beam, 100, 10, false);
%
%   % 指定缓存组名称（如收敛率测试）
%   wf_result = workflow.solve('PP', beam, 100, 10, false, 'convergence');

% 默认参数
if nargin < 5
    show_report = true;
end
if nargin < 6
    cache_prefix = '';
end

%% 输入验证
% 验证 beam 类型
assert(isa(beam, 'parameters.BeamModel'), ...
    'workflow:solve:InvalidBeamType', ...
    '参数 beam 必须是 parameters.BeamModel 对象，当前类型: %s', class(beam));

% 验证 BC
assert(ischar(BC) || isstring(BC), ...
    'workflow:solve:InvalidBCType', 'BC 必须是字符串');
BC = upper(char(BC));
assert(length(BC) == 2, ...
    'workflow:solve:InvalidBCLength', 'BC 必须是2个字符，如 PP, CF 等');
assert(ismember(BC, {'PP','CF','CP','PF','CC','FF','PR','CR','RR','RF'}), ...
    'workflow:solve:UnsupportedBC', ...
    '不支持的边界条件: %s\n支持的边界条件: PP, CF, CP, PF, CC, FF, PR, CR, RR, RF', BC);

% 验证 NElem
assert(isnumeric(NElem) && isscalar(NElem), ...
    'workflow:solve:InvalidNElemType', 'NElem 必须是数值标量');
assert(NElem > 0 && mod(NElem, 1) == 0, ...
    'workflow:solve:InvalidNElemValue', 'NElem 必须是正整数，当前值: %g', NElem);

% 验证 n_modes
assert(isnumeric(n_modes) && isscalar(n_modes), ...
    'workflow:solve:InvalidNModesType', 'n_modes 必须是数值标量');
assert(n_modes > 0 && mod(n_modes, 1) == 0, ...
    'workflow:solve:InvalidNModesValue', 'n_modes 必须是正整数，当前值: %g', n_modes);

% 内部日志函数：统一管理输出
    function log(varargin)
        if show_report
            fprintf(varargin{:});
        end
    end

% 物理形态列表
matrix_methods = {'axial', 'torsion', 'bending_shear'};
param_hash = beam.toHash();

%% 构建workflow缓存目录
% 缓存目录结构：
%   - 无前缀: .cache/{hash}_{BC}_N{NElem}_M{n_modes}/
%   - 有前缀: .cache/{cache_prefix}/{hash}_{BC}_N{NElem}_M{n_modes}/

% 获取项目根目录（src的上级目录）
script_path = mfilename('fullpath');
src_dir = fileparts(fileparts(script_path));  % 从 +workflow/solve.m 向上两级
project_root = fileparts(src_dir);  % src的上级目录

% 构建缓存名称
cache_name = sprintf('%s_%s_N%d_M%d', param_hash(1:8), BC, NElem, n_modes);

% 构建完整路径
if isempty(cache_prefix)
    workflow_root = fullfile(project_root, '.cache', cache_name);
else
    workflow_root = fullfile(project_root, '.cache', cache_prefix, cache_name);
end

if ~exist(workflow_root, 'dir')
    mkdir(workflow_root);
end

% 显示缓存位置
if isempty(cache_prefix)
    log('\n✓ 缓存目录: .cache/%s\n', cache_name);
else
    log('\n✓ 缓存目录: .cache/%s/%s\n', cache_prefix, cache_name);
end

%% 循环求解所有物理形态
matrix_methods = {'axial', 'torsion', 'bending_shear'};
for i = 1:length(matrix_methods)
    method = matrix_methods{i};

    % 创建子目录
    method_dir = fullfile(workflow_root, method);
    if ~exist(method_dir, 'dir')
        mkdir(method_dir);
    end

    fem_file = fullfile(method_dir, 'fem_result.mat');
    ana_file = fullfile(method_dir, 'ana_result.mat');

    % FEM求解（带缓存）
    if exist(fem_file, 'file')
        fem = load(fem_file);
        log('✓ 加载 FEM 缓存: %s\n', method);
    else
        fem = solvers.solve(BC, beam, method, NElem, n_modes);
        save(fem_file, '-struct', 'fem', '-v7.3');
        log('✓ 保存 FEM: %s\n', method);
    end

    % 解析解求解（带缓存）
    if exist(ana_file, 'file')
        ana = load(ana_file);
        log('✓ 加载解析解缓存: %s\n', method);
    else
        switch method
            case 'axial'
                ana = analytical.axial.solve(beam, BC, n_modes);
            case 'torsion'
                ana = analytical.torsion.solve(beam, BC, n_modes);
            case 'bending_shear'
                ana = analytical.bending_shear.solve(beam, BC, n_modes);
        end
        save(ana_file, '-struct', 'ana', '-v7.3');
        log('✓ 保存解析解: %s\n', method);
    end

    % 对比（自动打印）
    comp = comparison.compare(fem, ana, 'n_modes', n_modes, 'verbose', show_report);

    % 保存对比结果
    comp_file = fullfile(method_dir, 'comp_result.mat');
    save(comp_file, '-struct', 'comp', '-v7.3');
    log('✓ 保存对比: %s\n', method);

    % 打包结果
    wf_result.(method).fem = fem;
    wf_result.(method).ana = ana;
    wf_result.(method).comp = comp;
end

% 保存元数据
% 确保workflow_root目录存在（加载缓存时可能只有子目录）
if ~exist(workflow_root, 'dir')
    mkdir(workflow_root);
end

metadata.BC = BC;
metadata.NElem = NElem;
metadata.n_modes = n_modes;
metadata.beam_hash = param_hash;
metadata.beam_L = beam.L;
metadata.beam_E = beam.E;
metadata.beam_section_type = beam.section_type;
metadata.timestamp = datetime('now');
save(fullfile(workflow_root, 'metadata.mat'), '-struct', 'metadata', '-v7.3');

wf_result.save_path = workflow_root;

log('\n✓ Workflow 完成: %s\n', workflow_root);


%% 可选：输出报告
if show_report
    % 打印总体摘要
    printOverallSummary(BC, beam, NElem, n_modes, wf_result);
end

end

%% 辅助函数：打印总体摘要
function printOverallSummary(BC, beam, NElem, n_modes, wf_result)
% printOverallSummary - 打印总体参数摘要和误差汇总

fprintf('\n');
fprintf('========================================\n');
fprintf('WORKFLOW 求解完成\n');
fprintf('========================================\n\n');

% 输入参数摘要
fprintf('【求解配置】\n');
fprintf('  边界条件: %s\n', BC);
fprintf('  单元数:   %d\n', NElem);
fprintf('  模态数:   %d\n', n_modes);
fprintf('\n');

% 梁参数（按逻辑顺序：几何 → 截面 → 材料 → 计算属性）
fprintf('【梁参数】\n');
fprintf('  长度 L = %.3f m\n', beam.L);
fprintf('\n');

fprintf('  截面类型: %s\n', beam.section_type);
% 显示截面尺寸（根据类型动态显示）
if strcmp(beam.section_type, 'rectangular')
    fprintf('    宽度 b = %.3f m\n', beam.dims.b);
    fprintf('    高度 h = %.3f m\n', beam.dims.h);
elseif strcmp(beam.section_type, 'circular')
    fprintf('    直径 d = %.3f m\n', beam.dims.d);
elseif strcmp(beam.section_type, 'custom')
    fprintf('    自定义截面\n');
end
fprintf('\n');

fprintf('  材料属性:\n');
fprintf('    弹性模量 E   = %.2e Pa\n', beam.E);
fprintf('    剪切模量 G   = %.2e Pa\n', beam.G);
fprintf('    泊松比 ν     = %.3f\n', beam.nu);
fprintf('    密度 ρ       = %.1f kg/m³\n', beam.rho);
fprintf('\n');

fprintf('  截面属性:\n');
fprintf('    面积 A       = %.4e m²\n', beam.A);
fprintf('    惯性矩 I     = %.4e m⁴\n', beam.I);
fprintf('    极惯性矩 J   = %.4e m⁴\n', beam.J);
fprintf('    剪切系数 κ   = %.4f\n', beam.kappa);
fprintf('\n');

% 总体误差汇总
fprintf('【总体误差汇总】\n');
methods = {'axial', 'torsion', 'bending_shear'};
for i = 1:length(methods)
    method = methods{i};
    comp = wf_result.(method).comp;
    fprintf('  %-14s: 平均误差 %.4f%%, 最大误差 %.4f%%\n', ...
        upper(method), comp.summary.mean_error, comp.summary.max_error);
end
fprintf('\n');

end

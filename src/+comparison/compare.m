function comp = compare(fem_result, ana_result, varargin)
% compare - 对比有限元结果和解析解
%
% 用法:
%   comp = comparison.compare(fem_result, ana_result)
%   comp = comparison.compare(fem_result, ana_result, 'n_modes', 5)
%   comp = comparison.compare(fem_result, ana_result, 'verbose', true)
%
% 输入:
%   fem_result - 有限元求解结果 (solvers.solve 返回)
%   ana_result - 解析解结果 (analytical.*.solve 返回)
%
% 参数对:
%   'n_modes'  - 对比的模态数 (默认: min(fem, ana 模态数))
%   'verbose'  - 是否打印对比结果 (默认: true)
%
% 输出:
%   comp - 对比结构体
%     .freq_fem    - FEM 频率 [Hz] (n×1)
%     .freq_ana    - 解析频率 [Hz] (n×1)
%     .freq_error  - 频率相对误差 [%] (n×1)
%     .mac         - MAC 值 (n×1)
%                    axial/torsion: 标量模态的MAC
%                    bending_shear: 向量模态[w;θ]的MAC
%     .modes_fem   - FEM 模态（采样到节点）
%     .modes_ana   - 解析模态（采样到节点）
%     .summary     - 统计摘要
%       .mean_error - 平均频率误差 [%]
%       .max_error  - 最大频率误差 [%]
%       .mean_mac   - 平均 MAC
%       .min_mac    - 最小 MAC
%
% 示例:
%   % 弯曲剪切
%   fem_result = solvers.solve('PP', params, 'bending_shear');
%   ana_result = analytical.bending_shear.solve(params, 'PP', 10);
%   comp = comparison.compare(fem_result, ana_result);
%
%   % 轴向
%   fem_result = solvers.solve('CC', params, 'axial');
%   ana_result = analytical.axial.solve(params, 'CC', 10);
%   comp = comparison.compare(fem_result, ana_result, 'verbose', false);

% 解析参数
p = inputParser;
addParameter(p, 'n_modes', min(length(fem_result.freq), length(ana_result.freq_Hz)));
addParameter(p, 'verbose', true);
parse(p, varargin{:});
n_modes = p.Results.n_modes;
verbose = p.Results.verbose;

% 验证边界条件一致
if ~strcmp(fem_result.BC, ana_result.BC)
    error('边界条件不一致: fem=%s, ana=%s', fem_result.BC, ana_result.BC);
end

% 获取物理类型
matrix_method = fem_result.matrix_method;

% 1. 智能模态匹配
if strcmp(matrix_method, 'bending_shear')
    match_info = comparison.matchModes(fem_result, ana_result);

    n_fem = length(fem_result.freq);
    comp.freq_fem = fem_result.freq(1:n_fem);
    comp.freq_ana = nan(n_fem, 1);
    comp.match_status = match_info.match_status;
    comp.match_info = match_info;  % 保存完整信息供打印

    % 使用匹配索引填充
    for i = 1:n_fem
        if ~isnan(match_info.ana_indices(i))
            comp.freq_ana(i) = ana_result.freq_Hz(match_info.ana_indices(i));
        end
    end
    comp.freq_error = abs(comp.freq_fem - comp.freq_ana) ./ comp.freq_ana * 100;

    % 统计（仅匹配模态）
    matched = comp.match_status == "Matched";
    comp.summary.mean_error = mean(comp.freq_error(matched), 'omitnan');
    comp.summary.max_error = max(comp.freq_error(matched));
else
    % axial/torsion: 简单索引
    % 检查实际模态数，防止索引越界
    n_fem_actual = length(fem_result.freq);
    n_ana_actual = length(ana_result.freq_Hz);
    n_compare = min([n_modes, n_fem_actual, n_ana_actual]);

    if n_compare < n_modes
        warning('comparison:compare:insufficientModes', ...
            'FEM返回模态数不足(%d<%d)，原因：单元数NElem过少导致DOF不足。建议增加单元数。将对比前%d个模态。', ...
            n_fem_actual, n_modes, n_compare);
    end

    comp.freq_fem = fem_result.freq(1:n_compare);
    comp.freq_ana = ana_result.freq_Hz(1:n_compare);
    comp.freq_error = abs(comp.freq_fem - comp.freq_ana) ./ comp.freq_ana * 100;
    comp.match_status = repmat("Matched", n_compare, 1);

    % 统计
    comp.summary.mean_error = mean(comp.freq_error);
    comp.summary.max_error = max(comp.freq_error);
end

% 2. 对比模态（计算MAC）
if strcmp(matrix_method, 'bending_shear')
    % 采样解析模态到 FEM 节点
    NNodes = fem_result.NElem + 1;
    L = fem_result.params.L;
    x_fem = linspace(0, L, NNodes)';
    z_fem = x_fem / L;  % 无量纲坐标

    % 采样所有解析模态
    n_ana = length(ana_result.lambda);
    lambda_cell = analytical.bending_shear.packLambdaForModes(...
        ana_result.lambda, ana_result.sys.par.alpha);

    % 采样
    U_analytical = zeros(NNodes, n_ana);
    phi_analytical = zeros(NNodes, n_ana);

    mode_idx = 1;
    for region = 1:3
        n_modes_region = length(lambda_cell{region});
        if n_modes_region == 0
            continue;
        end

        % 获取该 region 的 chebfun（矩阵值）
        U_cheb = ana_result.modes{region}{1};
        phi_cheb = ana_result.modes{region}{2};

        if ~isempty(U_cheb)
            % 一次性采样该 region 的所有模态
            U_mat = U_cheb(z_fem);      % NNodes × n_modes_region
            phi_mat = phi_cheb(z_fem);

            % 逐个提取模态
            for i = 1:n_modes_region
                if mode_idx > n_ana
                    break;
                end

                U_analytical(:, mode_idx) = U_mat(:, i);
                phi_analytical(:, mode_idx) = phi_mat(:, i);

                mode_idx = mode_idx + 1;
            end
        else
            % 跳过空的 region
            mode_idx = mode_idx + n_modes_region;
        end
    end

    % 提取 FEM 模态的 w 和 theta
    total_DOF = size(fem_result.modes, 1);
    w_fem = fem_result.modes(1:2:total_DOF, 1:n_fem);
    theta_fem = fem_result.modes(2:2:total_DOF, 1:n_fem);

    % 构造完整向量 [w; θ]
    modes_fem = [w_fem; theta_fem];
    modes_ana_full = [U_analytical; phi_analytical];

    % 计算MAC（使用匹配索引）
    comp.mac = nan(n_fem, 1);
    for i = 1:n_fem
        if comp.match_status(i) == "Matched"
            ana_idx = match_info.ana_indices(i);
            comp.mac(i) = comparison.computeMAC(modes_fem(:,i), modes_ana_full(:,ana_idx));
        end
    end

    % 保存模态数据
    comp.modes_fem = modes_fem;
    comp.modes_ana = modes_ana_full;

    % MAC 统计（仅匹配模态）
    matched = comp.match_status == "Matched";
    comp.summary.mean_mac = mean(comp.mac(matched), 'omitnan');
    comp.summary.min_mac = min(comp.mac(matched));

else
    % Axial 或 Torsion：标量系统
    NNodes = fem_result.NElem + 1;
    L = fem_result.params.L;
    x_fem = linspace(0, L, NNodes)';

    % 采样解析模态到 FEM 节点
    modes_ana = ana_result.getModes(x_fem);
    modes_ana = modes_ana(:, 1:n_compare);

    % 提取 FEM 模态
    modes_fem = fem_result.modes(:, 1:n_compare);

    % 计算 MAC
    comp.mac = zeros(n_compare, 1);
    for i = 1:n_compare
        comp.mac(i) = comparison.computeMAC(modes_fem(:,i), modes_ana(:,i));
    end

    % 保存模态数据
    comp.modes_fem = modes_fem;
    comp.modes_ana = modes_ana;

    % MAC 统计
    comp.summary.mean_mac = mean(comp.mac);
    comp.summary.min_mac = min(comp.mac);
end

% 3. 打印结果（可选）
if verbose
    if strcmp(matrix_method, 'bending_shear')
        printComparisonResults(comp, matrix_method, n_fem);
    else
        printComparisonResults(comp, matrix_method, n_modes);
    end
end

end

%% 辅助函数
function printComparisonResults(comp, matrix_method, n_modes)
% 打印对比结果

fprintf('\n【对比结果 - %s】\n', matrix_method);

% 表头
fprintf('%-6s | %-12s | %-12s | %-10s | %-8s | %-30s\n', ...
        '模态', 'FEM(Hz)', '解析(Hz)', '误差(%)', 'MAC', '状态');
fprintf('%s\n', repmat('-', 1, 88));

% 打印每个模态
for i = 1:n_modes
    if comp.match_status(i) == "Matched"
        fprintf('%-6d | %12.4f | %12.4f | %10.4f | %8.4f | Matched\n', ...
            i, comp.freq_fem(i), comp.freq_ana(i), comp.freq_error(i), comp.mac(i));
    else
        fprintf('%-6d | %12.4f | %12s | %10s | %8s | %s\n', ...
            i, comp.freq_fem(i), 'N/A', 'N/A', 'N/A', comp.match_status(i));
    end
end

% 统计摘要
fprintf('\n【统计摘要】\n');
n_matched = sum(comp.match_status == "Matched");
n_missing = sum(comp.match_status ~= "Matched");

if n_matched > 0
    fprintf('  匹配模态：%d 个\n', n_matched);
    fprintf('  频率误差：平均 %.4f%%, 最大 %.4f%%\n', ...
            comp.summary.mean_error, comp.summary.max_error);
    fprintf('  MAC：     平均 %.4f, 最小 %.4f\n', ...
            comp.summary.mean_mac, comp.summary.min_mac);
end

if n_missing > 0
    fprintf('\n  ⚠️  缺失模态：%d 个 (数值奇异性: λ≈α=%.1f)\n', ...
            n_missing, comp.match_info.alpha);
    fprintf('      FEM索引：[%s]\n', num2str(find(comp.match_status ~= "Matched")'));
end

fprintf('========================================\n\n');
end

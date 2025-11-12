function match_info = matchModes(fem_result, ana_result)
% matchModes - 智能模态匹配（向量比对 + λ≈α缺失检测）
%
% 输入:
%   fem_result - FEM结果（包含 .freq, .omega, .params）
%   ana_result - 解析解结果（包含 .freq_Hz, .lambda, .sys.par.alpha/.T）
%
% 输出:
%   match_info.ana_indices    - 匹配的解析解索引数组 [1, NaN, 2, 3, ...]（NaN表示缺失）
%   match_info.match_status   - 状态数组 ["Matched", "Missing (Numerical Singularity)", ...]
%   match_info.alpha          - α值（供打印）
%
% 算法:
%   1. 贪婪匹配：从前往后逐个匹配
%   2. 匹配成功：ana_idx++
%   3. 匹配失败且λ≈α：标记缺失，ana_idx不变（下个FEM继续匹配当前解析解）
%   4. 匹配失败且λ不在α附近：报错

n_fem = length(fem_result.freq);
n_ana = length(ana_result.freq_Hz);

% 获取α和T
alpha = ana_result.sys.par.alpha;
T = ana_result.sys.par.T;

% 计算FEM的lambda
lambda_fem = (fem_result.omega .* T).^2;
lambda_ana = ana_result.lambda;

% 初始化
ana_indices = nan(n_fem, 1);
match_status = repmat("Matched", n_fem, 1);

% 全局匹配算法：对每个FEM模态，在ana中找最接近的λ
lambda_tolerance = 0.10;  % 10%相对误差容忍度（考虑高阶模态FEM误差增大）
alpha_tolerance = 0.10;   % 10%相对误差（判断λ≈α，仅标记非常接近α的模态）

used_ana_indices = false(n_ana, 1);  % 记录已使用的ana索引

for i = 1:n_fem
    % 优先检查：FEM模态是否接近α（解析解已知缺陷）
    is_near_alpha = abs(lambda_fem(i) - alpha) / alpha < alpha_tolerance;

    if is_near_alpha
        % FEM模态接近α：解析解在此处缺失
        ana_indices(i) = NaN;
        match_status(i) = "Missing (Numerical Singularity)";
    else
        % FEM模态不在α附近：尝试匹配解析解
        % 在ana_lambda中找最接近的λ（排除已使用的）
        relative_errors = abs(lambda_ana - lambda_fem(i)) ./ lambda_fem(i);
        relative_errors(used_ana_indices) = Inf;  % 排除已使用的ana模态
        [min_error, closest_idx] = min(relative_errors);

        if min_error < lambda_tolerance
            % 匹配成功：找到相近的λ
            ana_indices(i) = closest_idx;
            match_status(i) = "Matched";
            used_ana_indices(closest_idx) = true;
        else
            % 无法匹配：检查是否是FEM虚假模态
            max_ana_lambda = max(lambda_ana);

            if lambda_fem(i) > max_ana_lambda
                % FEM虚假模态：λ_fem超出解析解范围
                ana_indices(i) = NaN;
                match_status(i) = "Spurious (Increase NElem)";
            else
                % 真正的匹配失败 → 严重问题
                error('comparison:matchModes:unmatchedMode', ...
                      'FEM模态%d无法匹配：最小λ误差%.1f%%，λ_fem=%.2f, α=%.2f\n解析解范围: λ∈[%.2f, %.2f]', ...
                      i, min_error*100, lambda_fem(i), alpha, min(lambda_ana), max_ana_lambda);
            end
        end
    end
end

% 打包结果
match_info.ana_indices = ana_indices;
match_info.match_status = match_status;
match_info.alpha = alpha;

end

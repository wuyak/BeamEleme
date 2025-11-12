function [lambda, V_r, n_rigid] = filterEigen(lambda, V_r, BC, matrix_method, params, num_modes)
% filterEigen - 过滤和排序特征值
%
% 输入:
%   lambda        - 特征值 (未过滤)
%   V_r           - 特征向量
%   BC            - 边界条件
%   matrix_method - 物理性质 ('axial', 'torsion', 'bending_shear')
%   params        - BeamParameters对象
%   num_modes     - 需要保留的模态数
%
% 输出:
%   lambda  - 过滤后的特征值 (已排序, 已限制到num_modes)
%   V_r     - 对应的特征向量
%   n_rigid - 识别到的刚体模态数
%
% 说明:
%   - FF边界：保留刚体模态
%   - 其他边界：过滤刚体模态
%   - 自适应阈值识别刚体/弹性模态分界点
%   - 根据物理性质判断刚体模态数量

% 过滤负特征值
valid_idx = lambda > -1e-6;
lambda = lambda(valid_idx);
V_r = V_r(:, valid_idx);

if strcmp(BC, 'FF')
    % ========== FF边界：识别并过滤刚体模态 ==========
    [lambda_sorted, sort_idx] = sort(lambda);

    % 根据物理性质判断刚体模态数
    % 自适应刚体模态识别：找到特征值的跳变点
    % 刚体模态和弹性模态之间通常有至少3个数量级的差异
    expected_rigid = solvers.getRigidModeCount(matrix_method, BC);

    if length(lambda_sorted) > expected_rigid
        % 计算相邻特征值的比值
        ratio = lambda_sorted(expected_rigid+1) / max(lambda_sorted(expected_rigid), 1e-20);

        if ratio > 1e3
            % 存在明显跳变，前 expected_rigid 个是刚体模态
            rigid_threshold = sqrt(lambda_sorted(expected_rigid) * lambda_sorted(expected_rigid+1));
        else
            % 没有明显跳变，使用固定阈值
            rigid_threshold = 1e-3;
        end
    else
        rigid_threshold = 1e-3;
    end

    n_rigid = sum(lambda_sorted < rigid_threshold);

    % 过滤刚体模态（只保留弹性模态）
    elastic_idx = lambda_sorted >= rigid_threshold;
    lambda = lambda_sorted(elastic_idx);
    V_r = V_r(:, sort_idx(elastic_idx));

else
    % ========== 其他边界：过滤刚体模态 ==========
    % 基于最大特征值的相对大小设定自适应阈值
    lambda_max = max(abs(lambda));

    if lambda_max > 1e-6
        % 自适应阈值：相对于最大特征值的1e-10倍
        threshold = lambda_max * 1e-10;
    else
        % 所有特征值都很小，使用绝对阈值
        threshold = 1e-8;
    end

    valid_idx = lambda > threshold;
    lambda = lambda(valid_idx);
    V_r = V_r(:, valid_idx);

    n_rigid = 0;  % 已过滤
end

% 排序（FF边界已排序，其他边界需要排序）
if ~strcmp(BC, 'FF')
    [lambda, sort_idx] = sort(lambda);
    V_r = V_r(:, sort_idx);
end
% 限制模态数
n_modes = min(num_modes, length(lambda));
lambda = lambda(1:n_modes);
V_r = V_r(:, 1:n_modes);

end

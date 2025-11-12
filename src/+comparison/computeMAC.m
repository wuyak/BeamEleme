function mac_value = computeMAC(mode1, mode2)
% computeMAC - 计算两个模态的MAC值（Modal Assurance Criterion）
%
% MAC值衡量两个模态振型的相似度，范围0-1:
%   1.0  - 完全相同（理想）
%   >0.9 - 高度相关（很好）
%   0.7-0.9 - 中等相关（可接受）
%   <0.7 - 弱相关或不相关（差）
%
% 公式:
%   MAC(φ_i, φ_j) = |φ_i^T * φ_j|^2 / ((φ_i^T * φ_i) * (φ_j^T * φ_j))
%
% 用法:
%   mac = computeMAC(mode1, mode2)           % 两个列向量
%   mac_matrix = computeMAC(modes_A, modes_B) % 两个模态矩阵
%
% 输入:
%   mode1 - 第一个模态向量(列向量) 或 模态矩阵(DOF × n_modes)
%   mode2 - 第二个模态向量(列向量) 或 模态矩阵(DOF × n_modes)
%
% 输出:
%   mac_value - 标量MAC值 或 MAC矩阵 (n_modes1 × n_modes2)
%
% 示例:
%   % 单个模态对比
%   result1 = solveBeam('PP', 'nelems', 20);
%   result2 = solveBeam('PP', 'nelems', 100);
%   mac = computeMAC(result1.modes(:,1), result2.modes(:,1));
%   fprintf('MAC值: %.4f\n', mac);
%
%   % 多模态对比矩阵
%   mac_matrix = computeMAC(result1.modes(:,1:5), result2.modes(:,1:5));
%   imagesc(mac_matrix); colorbar;

    % 确保是列向量或矩阵
    if isrow(mode1)
        mode1 = mode1(:);
    end
    if isrow(mode2)
        mode2 = mode2(:);
    end

    % 检查维度
    if size(mode1, 1) ~= size(mode2, 1)
        error('两个模态的自由度数必须相同 (mode1: %d, mode2: %d)', ...
              size(mode1,1), size(mode2,1));
    end

    % 获取模态数
    n_modes1 = size(mode1, 2);
    n_modes2 = size(mode2, 2);

    % 计算MAC矩阵
    mac_value = zeros(n_modes1, n_modes2);

    for i = 1:n_modes1
        for j = 1:n_modes2
            phi_i = mode1(:, i);
            phi_j = mode2(:, j);

            % 计算内积
            numerator = abs(phi_i' * phi_j)^2;
            denominator = (phi_i' * phi_i) * (phi_j' * phi_j);

            % 避免除零
            if denominator < 1e-15
                mac_value(i, j) = 0;
            else
                mac_value(i, j) = numerator / denominator;
            end
        end
    end

    % 如果是单个模态，返回标量
    if n_modes1 == 1 && n_modes2 == 1
        mac_value = mac_value(1, 1);
    end
end

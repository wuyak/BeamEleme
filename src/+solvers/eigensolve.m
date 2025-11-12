function [lambda, V_r] = eigensolve(K_r, M_r, BC, matrix_method, n_modes)
% eigensolve - 求解广义特征值问题
%
% 输入:
%   K_r           - 缩减刚度矩阵
%   M_r           - 缩减质量矩阵
%   BC            - 边界条件 (用于判断是否病态)
%   matrix_method - 物理性质 ('axial', 'torsion', 'bending_shear')
%   n_modes       - 需要的模态数
%
% 输出:
%   lambda - 特征值 (未排序, 未过滤)
%   V_r    - 特征向量 (对应缩减自由度)
%
% 说明:
%   - 自动选择稀疏/密集求解器
%   - 根据物理性质和边界条件判断刚体模态数
%   - FF强制使用eig（多刚体模态，eigs数值不稳定）
%   - RF/RR在n_dofs>=150时使用eigs（可能有刚体模态，但数量少）
%   - eigs失败自动回退到eig

n_dofs = length(K_r);

% 根据物理性质和边界条件判断刚体模态数
% 使用统一的物理性质分类系统
n_rigid_expected = solvers.getRigidModeCount(matrix_method, BC);

% 识别病态边界条件（包含刚体模态）
% eigs对秩亏矩阵数值不稳定，有刚体模态时强制使用eig
is_ill_conditioned = (n_rigid_expected > 0);

% 决定使用密集求解器的条件
use_dense_solver = is_ill_conditioned || n_dofs < 100;

% 刚体模态安全余量
if is_ill_conditioned
    % 根据实际刚体模态数设置余量
    max_rigid = n_rigid_expected + 2;  % 预期刚体模态数 + 安全缓冲
else
    max_rigid = 2;  % 其他边界条件：保守余量
end

if ~use_dense_solver
    % 大规模问题：使用稀疏求解器
    % 策略：请求额外模态以覆盖可能的刚体模态

    % 请求的特征值数：目标模态数 + 安全余量
    k_eigs = min(n_modes + max_rigid, n_dofs - 1);

    % 暂时禁用 eigs 的病态矩阵警告
    warn_state = warning('query', 'MATLAB:eigs:IllConditionedA');
    warning('off', 'MATLAB:eigs:IllConditionedA');

    try
        [V_r, lambda_diag] = eigs(K_r, M_r, k_eigs, 'smallestabs');
        lambda = diag(lambda_diag);
    catch ME
        % eigs 失败，回退到密集求解器
        solvers.handleEigsFailure(ME, K_r, M_r);
        [V_r, lambda_diag] = eig(full(K_r), full(M_r));
        lambda = diag(lambda_diag);
    end

    % 恢复警告状态
    warning(warn_state.state, 'MATLAB:eigs:IllConditionedA');
else
    % 小规模问题或病态边界条件：使用密集求解器
    [V_r, lambda_diag] = eig(full(K_r), full(M_r));
    lambda = diag(lambda_diag);
end

end

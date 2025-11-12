function handleEigsFailure(ME, K, M)
% handleEigsFailure - 处理eigs求解器失败
%
% 输入:
%   ME - MATLAB异常对象
%   K  - 刚度矩阵
%   M  - 质量矩阵

% 记录矩阵信息以便调试
matrix_size = size(K);
cond_K = cond(full(K));
cond_M = cond(full(M));

warning('solvers:eigs_failed', ...
        'eigs求解器失败: %s\n矩阵信息: %dx%d, cond(K)=%.2e, cond(M)=%.2e\n回退到密集eig求解器', ...
        ME.message, matrix_size(1), cond_K, cond_M);

end
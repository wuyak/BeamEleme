function freq_Hz = getTimoshenkoAnalyticalFrequencies(params, BC, num_modes)
% getTimoshenkoAnalyticalFrequencies - 获取Timoshenko梁的解析频率
%
% 输入:
%   params - BeamParameters对象，包含梁的几何和材料参数
%   BC - 边界条件字符串:
%        "PP": 简支-简支 (pinned-pinned)
%        "CF": 固支-自由 (clamped-free, 悬臂)
%        "CC": 固支-固支 (clamped-clamped)
%        "FF": 自由-自由 (free-free)
%        等等 (见getTimoshenkoEigenvalues.m)
%   num_modes - 要计算的模态数
%
% 输出:
%   freq_Hz - 解析频率 (Hz), 列向量
%
% 理论来源:
%   Khasawneh, F. A., & Segalman, D. (2019).
%   Exact and numerically stable expressions for Euler-Bernoulli
%   and Timoshenko beam modes. Applied Acoustics, 151, 215-228.
%
% 依赖:
%   - chebfun工具箱
%   - getTimoshenkoEigenvalues.m

% 添加chebfun路径（如果还没有添加）
chebfun_path = fullfile(fileparts(mfilename('fullpath')), '..', '..', '..', 'chebfun');
if exist(chebfun_path, 'dir') && ~exist('chebfun', 'file')
    addpath(chebfun_path);
    chebfunroot;
end

% 构建系统参数
sys = analytical.bending_shear.buildSystemParameters(params);

% 获取特征值的chebfun表示
eigVal_chebfun = analytical.bending_shear.getEigenvalues(sys, BC);

% 求解lambda
lambda_vals = analytical.bending_shear.solveLambdaFromChebfun(eigVal_chebfun, sys.par.alpha);

if length(lambda_vals) < num_modes
    warning('只找到 %d 个模态，少于请求的 %d 个', length(lambda_vals), num_modes);
    num_modes = length(lambda_vals);
end

lambda_vals = lambda_vals(1:num_modes);

% 将lambda转换为圆频率omega (rad/s)
% 根据Khasawneh论文: lambda = (omega * T)^2
% 所以: omega = sqrt(lambda) / T
omega_rad_s = sqrt(lambda_vals) / sys.par.T;

% 转换为频率 (Hz)
freq_Hz = omega_rad_s / (2 * pi);

end

function freq_Hz = getTorsionAnalyticalFrequencies(params, BC, num_modes)
% getTorsionAnalyticalFrequencies - 获取扭转振动的解析频率
%
% 输入:
%   params - BeamParameters对象，包含梁的几何和材料参数
%   BC - 边界条件字符串:
%        "PP", "CF", "CC", "FF", "CP", "PF", "PR", "CR", "RR", "RF"
%   num_modes - 要计算的模态数
%
% 输出:
%   freq_Hz - 解析频率 (Hz), 列向量
%
% 理论来源:
%   扭转振动通用公式:
%   f_i = (λ_i / 2πL) * sqrt(JG / ρIp)
%
%   边界条件映射（[u,w,φ,θ] -> Torsion [φ]）:
%   - C, R -> C (约束φ)
%   - P, F -> F (不约束φ)
%
% 特征值 λ_i:
%   - CC: i*π       (i = 1,2,3,...)
%   - CF: (2i-1)π/2 (i = 1,2,3,...)
%   - FC: (2i-1)π/2 (i = 1,2,3,...)
%   - FF: i*π       (i = 1,2,3,...)

% 验证边界条件
assert(length(BC) == 2, '边界条件必须是2字符，如PP, CF等');
BC = upper(BC);

% 从BeamParameters提取物理参数
G = params.G;           % 剪切模量 (Pa)
rho = params.rho;       % 密度 (kg/m³)
L = params.L;           % 梁长 (m)
J = params.J;           % 扭转常数 (m⁴)
Ip = params.Ip;         % 极面积惯性矩 (m⁴)

% 当前2d梁，扭转的边界条件只有 C 和 F 两种，需要建立CPRF到CF的映射
torsion_BC = analytical.torsion.mapBoundaryConditionToTorsion(BC);

% 根据边界条件计算特征值 λ_i
switch torsion_BC
    case "CC"  % 固定-固定
        % λ_i = i*π
        i_values = (1:num_modes)';
        lambda_i = i_values * pi;

    case "CF"  % 固定-自由
        % λ_i = (2i-1)π/2
        i_values = (1:num_modes)';
        lambda_i = (2*i_values - 1) * pi / 2;

    case "FC"  % 自由-固定
        % λ_i = (2i-1)π/2
        i_values = (1:num_modes)';
        lambda_i = (2*i_values - 1) * pi / 2;

    case "FF"  % 自由-自由
        % λ_i = i*π (去除刚体模态 i=0)
        i_values = (1:num_modes)';
        lambda_i = i_values * pi;

    otherwise
        error('不支持的扭转边界条件: %s', torsion_BC);
end

% 计算频率 (Hz)
% f_i = (λ_i / 2πL) * sqrt(JG / ρIp)
freq_Hz = (lambda_i / (2*pi*L)) * sqrt(J*G / (rho*Ip));

end

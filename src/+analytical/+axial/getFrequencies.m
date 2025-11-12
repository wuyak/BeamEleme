function freq_Hz = getAxialAnalyticalFrequencies(params, BC, num_modes)
% getAxialAnalyticalFrequencies - 获取轴向（压缩）振动的解析频率
%
% 输入:
%   params - BeamParameters对象，包含梁的几何和材料参数
%   BC - 边界条件字符串（4DOF定义）:
%        "PP", "CF", "CC", "FF", "CP", "PF", "PR", "CR", "RR", "RF"
%   num_modes - 要计算的模态数
%
% 输出:
%   freq_Hz - 解析频率 (Hz), 列向量
%
% 理论来源:
%   压缩振动通用公式:
%   f_i = (λ_i / 2πL) * sqrt(E / ρ)
%
%   边界条件映射（ [u,w,φ,θ] -> Axial [u]）:
%   - C, P, R -> C (约束u)
%   - F -> F (不约束u)
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
E = params.E;           % 杨氏模量 (Pa)
rho = params.rho;       % 密度 (kg/m³)
L = params.L;           % 梁长 (m)

% 当前2d梁，轴向的边界条件只有 C 和 F 两种，需要建立CPRF到CF的映射
axial_BC = analytical.axial.mapBoundaryConditionToAxial(BC);

% 根据边界条件计算特征值 λ_i
switch axial_BC
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
        error('不支持的轴向边界条件: %s', axial_BC);
end

% 计算频率 (Hz)
% f_i = (λ_i / 2πL) * sqrt(E / ρ)
freq_Hz = (lambda_i / (2*pi*L)) * sqrt(E / rho);

end

function phi_modes = getTorsionModes(x, L, BC, num_modes)
% getTorsionModes - 获取扭转振动的解析模态形状
%
% 输入:
%   x - 空间坐标 (m)，列向量
%   L - 梁长 (m)
%   BC - 边界条件字符串（4DOF定义）:
%        "PP", "CF", "CC", "FF", "CP", "PF", "PR", "CR", "RR", "RF"
%   num_modes - 要计算的模态数
%
% 输出:
%   phi_modes - 扭转角模态 [N_nodes × num_modes]
%
% 理论来源:
%   扭转振动模态形状（归一化）:
%   - CC: φ_i(x) = sin(iπx/L)
%   - CF: φ_i(x) = sin((2i-1)πx/(2L))
%   - FC: φ_i(x) = cos((2i-1)πx/(2L))
%   - FF: φ_i(x) = cos(iπx/L)

% 验证边界条件
assert(length(BC) == 2, '边界条件必须是2字符，如PP, CF等');
BC = upper(BC);

% 映射4DOF边界条件到扭转1D边界条件
torsion_BC = analytical.torsion.mapBoundaryConditionToTorsion(BC);

% 初始化模态矩阵
phi_modes = zeros(length(x), num_modes);

% 根据边界条件计算模态形状
switch torsion_BC
    case "CC"  % 固定-固定
        for i = 1:num_modes
            phi_modes(:, i) = sin(i * pi * x / L);
        end

    case "CF"  % 固定-自由
        for i = 1:num_modes
            phi_modes(:, i) = sin((2*i - 1) * pi * x / (2*L));
        end

    case "FC"  % 自由-固定
        for i = 1:num_modes
            phi_modes(:, i) = cos((2*i - 1) * pi * x / (2*L));
        end

    case "FF"  % 自由-自由
        for i = 1:num_modes
            phi_modes(:, i) = cos(i * pi * x / L);
        end

    otherwise
        error('不支持的扭转边界条件: %s', torsion_BC);
end

% 归一化（单位模）
for i = 1:num_modes
    phi_modes(:, i) = phi_modes(:, i) / norm(phi_modes(:, i));
end

end

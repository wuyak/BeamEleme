function u_modes = getAxialModes(x, L, BC, num_modes)
% getAxialModes - 获取轴向振动的解析模态形状
%
% 输入:
%   x - 空间坐标 (m)，列向量
%   L - 梁长 (m)
%   BC - 边界条件字符串（4DOF定义）:
%        "PP", "CF", "CC", "FF", "CP", "PF", "PR", "CR", "RR", "RF"
%   num_modes - 要计算的模态数
%
% 输出:
%   u_modes - 轴向位移模态 [N_nodes × num_modes]
%
% 理论来源:
%   轴向振动模态形状（归一化）:
%   - CC: u_i(x) = sin(iπx/L)
%   - CF: u_i(x) = sin((2i-1)πx/(2L))
%   - FC: u_i(x) = cos((2i-1)πx/(2L))
%   - FF: u_i(x) = cos(iπx/L)

% 验证边界条件
assert(length(BC) == 2, '边界条件必须是2字符，如PP, CF等');
BC = upper(BC);

% 映射4DOF边界条件到轴向1D边界条件
axial_BC = analytical.axial.mapBoundaryConditionToAxial(BC);

% 初始化模态矩阵
u_modes = zeros(length(x), num_modes);

% 根据边界条件计算模态形状
switch axial_BC
    case "CC"  % 固定-固定
        for i = 1:num_modes
            u_modes(:, i) = sin(i * pi * x / L);
        end

    case "CF"  % 固定-自由
        for i = 1:num_modes
            u_modes(:, i) = sin((2*i - 1) * pi * x / (2*L));
        end

    case "FC"  % 自由-固定
        for i = 1:num_modes
            u_modes(:, i) = cos((2*i - 1) * pi * x / (2*L));
        end

    case "FF"  % 自由-自由
        for i = 1:num_modes
            u_modes(:, i) = cos(i * pi * x / L);
        end

    otherwise
        error('不支持的轴向边界条件: %s', axial_BC);
end

% 归一化（单位模）
for i = 1:num_modes
    u_modes(:, i) = u_modes(:, i) / norm(u_modes(:, i));
end

end

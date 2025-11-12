function axial_BC = mapBoundaryConditionToAxial(BC_4DOF)
% mapBoundaryConditionToAxial - 将4DOF边界条件映射到轴向1D边界条件
%
% 输入:
%   BC_4DOF - 边界条件字符串，如 "PP", "CF", "CC" 等
%
% 输出:
%   axial_BC - 轴向边界条件 ("CC", "CF", "FF")
%
% 边界条件映射规则:
%   边界条件 [u,w,φ,θ] -> 轴向 [u]:
%   - C, P, R -> C (约束u)
%   - F -> F (不约束u)
%
% 示例:
%   axial_BC = mapBoundaryConditionToAxial("CF")  % 返回 "CF"
%   axial_BC = mapBoundaryConditionToAxial("PR")  % 返回 "CC"

% 验证输入
assert(length(BC_4DOF) == 2, '边界条件必须是2字符，如PP, CF等');
BC_4DOF = upper(BC_4DOF);

% 提取左右端边界条件
left = BC_4DOF(1);
right = BC_4DOF(2);

% 映射函数：C,P,R -> C; F -> F
map_axial = @(bc) char(ismember(bc, ['C', 'P', 'R']) * 'C' + ~ismember(bc, ['C', 'P', 'R']) * 'F');

% 应用映射
axial_BC = [map_axial(left), map_axial(right)];

% 验证是否为支持的边界条件
supported_BCs = ["CC", "CF", "FC", "FF"];
if ~ismember(axial_BC, supported_BCs)
    error('轴向边界条件 %s 未实现。支持的边界条件: %s', ...
        axial_BC, strjoin(supported_BCs, ', '));
end

end
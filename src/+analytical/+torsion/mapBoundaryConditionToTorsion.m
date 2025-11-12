function torsion_BC = mapBoundaryConditionToTorsion(BC_4DOF)
% mapBoundaryConditionToTorsion - 将4DOF边界条件映射到扭转1D边界条件
%
% 输入:
%   BC_4DOF - 4DOF边界条件字符串，如 "PP", "CF", "CC" 等
%
% 输出:
%   torsion_BC - 扭转1D边界条件 ("CC", "CF", "FF")
%
% 边界条件映射规则:
%   4DOF边界条件 [u,w,φ,θ] -> 扭转 [φ]:
%   - C, R -> C (约束φ)
%   - P, F -> F (不约束φ)
%
% 示例:
%   torsion_BC = mapBoundaryConditionToTorsion("CF")  % 返回 "CF"
%   torsion_BC = mapBoundaryConditionToTorsion("PR")  % 返回 "CF"

% 验证输入
assert(length(BC_4DOF) == 2, '边界条件必须是2字符，如PP, CF等');
BC_4DOF = upper(BC_4DOF);

% 提取左右端边界条件
left = BC_4DOF(1);
right = BC_4DOF(2);

% 映射函数：C,R -> C; P,F -> F
map_torsion = @(bc) char(ismember(bc, ['C', 'R']) * 'C' + ~ismember(bc, ['C', 'R']) * 'F');

% 应用映射
torsion_BC = [map_torsion(left), map_torsion(right)];

% 验证是否为支持的边界条件
supported_BCs = ["CC", "CF", "FC", "FF"];
if ~ismember(torsion_BC, supported_BCs)
    error('扭转边界条件 %s 未实现。支持的边界条件: %s', ...
        torsion_BC, strjoin(supported_BCs, ', '));
end

end
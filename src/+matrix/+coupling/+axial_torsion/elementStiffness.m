"function Ke = elementStiffness(params)
% elementStiffness - 轴向-扭转耦合单元刚度矩阵（已弃用）
%
% 警告: 这是混合求解方案，已被轴向/扭转分离求解替代
% 建议使用: matrix.stiffness.axial.element 和 matrix.stiffness.torsion.element
%
% 输入:
%   params - 结构体，包含材料、几何参数
%
% 输出:
%   Ke - 4x4 单元刚度矩阵
%        DOFs: [u1, phi1, u2, phi2]
%        u:   轴向位移
%        phi: 扭转角

% 提取参数
L = params.dL;  % 单元长度
E = params.E;
G = params.G;
A = params.A;
J = params.J;

% 轴向刚度系数
k_axial = E * A / L;

% 扭转刚度系数
k_torsion = G * J / L;

% 组装块对角矩阵（轴向与扭转解耦）
Ke = zeros(4, 4);

% 轴向部分 (DOF 1, 3)
Ke([1, 3], [1, 3]) = k_axial * [1, -1; -1, 1];

% 扭转部分 (DOF 2, 4)
Ke([2, 4], [2, 4]) = k_torsion * [1, -1; -1, 1];

end"
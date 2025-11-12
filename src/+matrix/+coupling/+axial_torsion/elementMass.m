"function Me = elementMass(params)
% elementMass - 轴向-扭转耦合单元质量矩阵（已弃用）
%
% 警告: 这是混合求解方案，已被轴向/扭转分离求解替代
% 建议使用: matrix.mass.axial.element 和 matrix.mass.torsion.element
%
% 输入:
%   params - 结构体，包含材料、几何参数
%
% 输出:
%   Me - 4x4 单元质量矩阵
%        DOFs: [u1, phi1, u2, phi2]
%        u:   轴向位移
%        phi: 扭转角

% 提取参数
L = params.dL;  % 单元长度
rho = params.rho;
A = params.A;
Ip = params.Ip;

% 轴向质量系数
m_axial = rho * A * L / 6;

% 扭转质量系数
m_torsion = rho * Ip * L / 6;

% 基础2x2一致质量矩阵（线性形函数）
M_consistent = [2, 1; 1, 2];

% 组装块对角矩阵
Me = zeros(4, 4);

% 轴向部分 (DOF 1, 3)
Me([1, 3], [1, 3]) = m_axial * M_consistent;

% 扭转部分 (DOF 2, 4)
Me([2, 4], [2, 4]) = m_torsion * M_consistent;

end"
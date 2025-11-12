"function [K, M] = assemble(params)
% assemble - 轴向-扭转耦合装配（已弃用）
%
% 警告: 这是混合求解方案，已被轴向/扭转分离求解替代
% 建议使用: solvers.assembly.assembleGlobalMatrices(params, 'axial')
%       solvers.assembly.assembleGlobalMatrices(params, 'torsion')
%
% 输入：
%   params - BeamParameters对象
%            必须满足: node_type = 2, DOF_per_node = 2
%
% 输出：
%   K - 全局刚度矩阵（稀疏）
%   M - 全局质量矩阵（稀疏）

% 输入验证
if params.node_type ~= 2
    error('轴向-扭转耦合模块仅支持2节点单元');
end

if params.DOF_per_node ~= 2
    error('轴向-扭转耦合模块: DOF_per_node必须为2 (u, phi)');
end

% 创建单元参数结构体
elem_params = struct();
elem_params.E = params.E;
elem_params.G = params.G;
elem_params.rho = params.rho;
elem_params.A = params.A;
elem_params.J = params.J;
elem_params.Ip = params.Ip;
elem_params.dL = params.L / params.NElem;

% 总自由度数
total_DOF = params.NNodes * params.DOF_per_node;
K = zeros(total_DOF, total_DOF);
M = zeros(total_DOF, total_DOF);

% 装配循环
for elem = 1:params.NElem
    % 计算单元矩阵
    Ke = matrix.coupling.axial_torsion.elementStiffness(elem_params);
    Me = matrix.coupling.axial_torsion.elementMass(elem_params);

    % 计算节点自由度索引
    node1 = elem;
    node2 = elem + 1;
    node1_dofs = (node1 - 1) * params.DOF_per_node + (1:params.DOF_per_node);
    node2_dofs = (node2 - 1) * params.DOF_per_node + (1:params.DOF_per_node);
    elem_dofs = [node1_dofs, node2_dofs];

    % 装配到全局矩阵
    K(elem_dofs, elem_dofs) = K(elem_dofs, elem_dofs) + Ke;
    M(elem_dofs, elem_dofs) = M(elem_dofs, elem_dofs) + Me;
end

end"
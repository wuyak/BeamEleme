function result = buildResult(lambda, V_r, free_dofs, BC, params, matrix_method, NElem)
% buildResult - 构建结果结构体
%
% 输入:
%   lambda        - 特征值
%   V_r           - 特征向量 (缩减自由度)
%   free_dofs     - 自由度索引
%   BC            - 边界条件
%   params        - BeamParameters对象
%   matrix_method - 矩阵装配方法
%   NElem         - 单元数
%
% 输出:
%   result - 结果结构体 (6个字段)
%     .freq          - 固有频率 [Hz]
%     .omega         - 角频率 [rad/s]
%     .modes         - 模态矩阵，物理意义取决于 matrix_method:
%                      Axial:        (N_nodes × n_modes) - 轴向位移 u
%                      Torsion:      (N_nodes × n_modes) - 扭转角 φ
%                      BendingShear: (2*N_nodes × n_modes) - 奇数行: w, 偶数行: θ
%     .BC            - 边界条件
%     .params        - 参数对象 (引用)
%     .matrix_method - 物理性质 ('axial' | 'torsion' | 'bending_shear')

% 计算频率
omega = sqrt(lambda);
freq = omega / (2*pi);

% 从 matrix_method 推断 DOF_per_node
DOF_per_node = solvers.getDOFPerNode(matrix_method);

% 重构完整模态向量
total_DOF = (NElem + 1) * DOF_per_node;
n_modes = size(V_r, 2);
modes = zeros(total_DOF, n_modes);
modes(free_dofs, :) = V_r;

% 封装结构体
result.freq = freq;
result.omega = omega;
result.modes = modes;
result.BC = BC;
result.params = params;
result.matrix_method = matrix_method;  % 保留此字段用于验证
result.NElem = NElem;  % 添加单元数字段

end

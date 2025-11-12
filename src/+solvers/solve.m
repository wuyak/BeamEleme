function fem_result = solve(BC, params, matrix_method, NElem, varargin)
% solve - 求解梁的固有频率和模态（有限元数值解）
%
% 用法:
%   fem_result = solvers.solve('PP', params, 'bending_shear', 100);
%   fem_result = solvers.solve('PP', params, 'bending_shear', 100, 10);
%   fem_result = solvers.solve('CC', params, 'axial', 50, 20);
%
% 输入:
%   BC     - 边界条件字符串
%            'PP' = Pinned-Pinned (简支-简支)
%            'CC' = Clamped-Clamped (固支-固支)
%            'CF' = Clamped-Free (悬臂)
%            'FF' = Free-Free (自由-自由)
%            'CP' = Clamped-Pinned
%            'PF' = Pinned-Free
%            'RR' = Roller-Roller (滚支-滚支)
%            'CR' = Clamped-Roller
%            'PR' = Pinned-Roller
%            'RF' = Roller-Free
%
%   params - BeamParameters对象
%            包含材料、几何参数
%
%   matrix_method - 物理性质类型
%            'axial'         - 轴向振动 (u: 轴向位移)
%            'torsion'       - 扭转振动 (φ: 扭转角)
%            'bending_shear' - 弯曲剪切振动 (w: 横向位移, θ: 转角)
%
%   NElem - 单元数
%
%   num_modes - (可选) 求解模态数，默认 10
%
% 输出:
%   fem_result - 数值解结构体 (6个字段)
%     .freq          - 固有频率 [Hz] (向量)
%     .omega         - 角频率 [rad/s] (向量)
%     .modes         - 模态矩阵（离散节点值），物理意义取决于 matrix_method:
%                      Axial:        (N_nodes × n_modes) - 轴向位移 u
%                      Torsion:      (N_nodes × n_modes) - 扭转角 φ
%                      BendingShear: (2*N_nodes × n_modes) - 奇数行: w, 偶数行: θ
%     .BC            - 边界条件 (字符串)
%     .params        - BeamParameters对象 (引用)
%     .matrix_method - 物理性质 ('axial' | 'torsion' | 'bending_shear')

% 示例:
%   % 分离求解模式
%   mat = parameters.MaterialLibrary.Steel();
%   sec = parameters.SectionLibrary.Rectangular(0.2, 0.1);
%   params = parameters.BeamParameters(mat, sec, 'L', 1.0);
%
%   % 轴向振动（数值解）
%   fem_result_axial = solvers.solve('CC', params, 'axial', 50);
%
%   % 扭转振动（数值解）
%   fem_result_torsion = solvers.solve('CC', params, 'torsion', 50, 20);
%
%   % 弯曲剪切振动（数值解）
%   fem_result_bending = solvers.solve('PP', params, 'bending_shear', 100);
%
% 说明:
%   - 分离模式: 每个求解只处理单一物理过程
%   - NElem 作为外部参数传入，支持收敛性分析
%   - num_modes 作为可选参数，默认 10

% 解析可选参数 num_modes
if nargin >= 5 && ~isempty(varargin{1})
    num_modes = varargin{1};
else
    num_modes = 10;  % 默认值
end

% 验证参数
params.validateParameters();

% 1. 装配全局矩阵
[K, M] = solvers.assemble(params, matrix_method, NElem);

% 2. 应用边界条件
[K_r, M_r, free_dofs] = boundary_conditions.BoundaryConditionManager.apply(...
    K, M, BC, params, matrix_method, NElem);

% 3. 求解特征系统
[lambda, V_r] = solvers.eigensolve(...
    K_r, M_r, BC, matrix_method, num_modes);

% 4. 过滤刚体模态
[lambda, V_r, ~] = solvers.filterEigen(...
    lambda, V_r, BC, matrix_method, params, num_modes);

% 5. 构建结果
fem_result = solvers.buildResult(...
    lambda, V_r, free_dofs, BC, params, matrix_method, NElem);

end
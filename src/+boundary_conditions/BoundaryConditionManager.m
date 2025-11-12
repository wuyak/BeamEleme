classdef BoundaryConditionManager
% BoundaryConditionManager - 边界条件管理器
%
% 功能：
%   - 支持分离模式的边界条件定义和应用
%   - 支持1-DOF(axial/torsion)和2-DOF(bending_shear)
%   - 根据物理量（u/φ/w/θ）判断边界条件约束
%
% 边界条件定义（基于4-DOF [u, w, φ, θ]）:
%   - C (Clamped):  约束 [u, w, φ, θ]
%   - P (Pinned):   约束 [u, w]
%   - R (Roller):   约束 [u, φ, θ]
%   - F (Free):     无约束
%
% 对应到1-DOF模式:
%   - axial模式 (u):   C/P/R → C; F → F
%   - torsion模式 (φ): C/R → C; P/F → F
% 对应到2-DOF模式:
%   - bending_shear模式 (w, θ):
%     C: w=0, θ=0
%     P: w=0, θ自由
%     R: w自由, θ自由
%     F: w自由, θ自由
%
% 用法：
%   [K_r, M_r, free_dofs] = BoundaryConditionManager.apply(K, M, BC, params, matrix_method);

    methods (Static)

        function [K_reduced, M_reduced, free_dofs] = apply(K, M, BC_String, params, matrix_method, NElem)
        % apply - 应用边界条件，返回缩减矩阵
        %
        % 输入:
        %   K, M          - 全局刚度和质量矩阵
        %   BC_String     - 边界条件字符串 (如'PP', 'CF')
        %   params        - BeamParameters对象
        %   matrix_method - 矩阵方法 ('axial', 'torsion', 'bending_shear')
        %   NElem         - 单元数
        %
        % 输出:
        %   K_reduced - 缩减后刚度矩阵
        %   M_reduced - 缩减后质量矩阵
        %   free_dofs - 自由度索引向量

            assert(length(BC_String) == 2, '边界条件必须是2字符，如PP, CF等');
            BC_String = upper(BC_String);
            assert(ismember(BC_String, {'PP','CF','CP','PF','CC','FF','PR','CR','RR','RF'}), ...
                '无效边界条件: %s\n支持的10种组合: PP,CF,CP,PF,CC,FF,PR,CR,RR,RF', BC_String);

            % 从 matrix_method 推断 DOF_per_node
            DOF_per_node = solvers.getDOFPerNode(matrix_method);

            % 获取约束自由度
            left_dofs = boundary_conditions.BoundaryConditionManager.getConstrainedDOFs(...
                BC_String(1), 1, DOF_per_node, matrix_method);
            right_dofs = boundary_conditions.BoundaryConditionManager.getConstrainedDOFs(...
                BC_String(2), NElem + 1, DOF_per_node, matrix_method);
            constrained_dofs = [left_dofs, right_dofs];

            % 计算自由度
            total_DOF = (NElem + 1) * DOF_per_node;
            free_dofs = setdiff(1:total_DOF, constrained_dofs);

            assert(~isempty(free_dofs), '过度约束: 所有自由度均被约束');

            % 矩阵缩减
            K_reduced = K(free_dofs, free_dofs);
            M_reduced = M(free_dofs, free_dofs);
        end

        function constrained_dofs = getConstrainedDOFs(BC_type, node, DOF_per_node, matrix_method)
        % getConstrainedDOFs - 获取边界条件约束的自由度（数据驱动）
        %
        % 真正按物理性质分类：使用映射表消除所有switch-case和重复代码
        %
        % 输入:
        %   BC_type       - 边界条件类型 ('C', 'P', 'F', 'R')
        %   node          - 节点编号
        %   DOF_per_node  - 每节点自由度数
        %   matrix_method - 物理性质 ('axial', 'torsion', 'bending_shear')
        %
        % 输出:
        %   constrained_dofs - 被约束的自由度索引向量

            base_dof = (node - 1) * DOF_per_node;

            % 定义物理量到约束规则的映射表
            % true = 被约束, false = 自由
            constraint_map = struct(...
                'u',     struct('C',true, 'P',true, 'R',true,  'F',false), ...
                'phi',   struct('C',true, 'P',false,'R',true,  'F',false), ...
                'w',     struct('C',true, 'P',true, 'R',false, 'F',false), ...
                'theta', struct('C',true, 'P',false,'R',true,  'F',false) ...
            );

            % 定义模式使用的物理量及其DOF位置
            mode_dofs = struct(...
                'axial',         struct('quantities', {{'u'}},            'positions', [1]), ...
                'torsion',       struct('quantities', {{'phi'}},          'positions', [1]), ...
                'bending_shear', struct('quantities', {{'w', 'theta'}},   'positions', [1, 2]) ...
            );

            % 获取当前模式的物理量配置
            if ~isfield(mode_dofs, matrix_method)
                error('未知的matrix_method: %s (应为 axial, torsion 或 bending_shear)', matrix_method);
            end
            config = mode_dofs.(matrix_method);
            physical_quantities = config.quantities;
            dof_positions = config.positions;

            % 根据物理量查表决定约束
            constrained_dofs = [];
            for i = 1:length(physical_quantities)
                qty = physical_quantities{i};
                if constraint_map.(qty).(BC_type)
                    constrained_dofs(end+1) = base_dof + dof_positions(i);
                end
            end
        end

    end
end

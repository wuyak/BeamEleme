function [K, M] = assemble(params, matrix_method, NElem)
% assemble - 装配全局刚度和质量矩阵
%
% 输入:
%   params        - BeamParameters对象
%   matrix_method - 'axial' | 'torsion' | 'bending_shear'
%   NElem         - 单元数
%
% 输出:
%   K - 全局刚度矩阵
%   M - 全局质量矩阵
%
% 说明:
%   - 按物理性质组装，物理矩阵已分类存储
%   - 仅支持分离模式，无耦合接口

    % 选择单元矩阵函数并确定DOF
    switch matrix_method
        case 'axial'
            stiffness_func = @matrix.stiffness.axial.element;
            mass_func = @matrix.mass.axial.element;
            DOF_per_node = 1;

        case 'torsion'
            stiffness_func = @matrix.stiffness.torsion.element;
            mass_func = @matrix.mass.torsion.element;
            DOF_per_node = 1;

        case 'bending_shear'
            stiffness_func = @matrix.stiffness.bending_shear.element;
            mass_func = @matrix.mass.bending_shear.element;
            DOF_per_node = 2;

        otherwise
            error('未知的 matrix_method: %s。仅支持 axial, torsion, bending_shear', matrix_method);
    end

    % 总自由度数
    total_DOF = (NElem + 1) * DOF_per_node;
    K = zeros(total_DOF, total_DOF);
    M = zeros(total_DOF, total_DOF);

    % 创建单元参数结构体 - 按需传递
    elem_params = struct();
    elem_params.L = params.L / NElem;

    switch matrix_method
        case 'axial'
            elem_params.E = params.E;
            elem_params.A = params.A;
            elem_params.rho = params.rho;

        case 'torsion'
            elem_params.G = params.G;
            elem_params.J = params.J;
            elem_params.Ip = params.Ip;
            elem_params.rho = params.rho;

        case 'bending_shear'
            elem_params.E = params.E;
            elem_params.I = params.I;
            elem_params.G = params.G;
            elem_params.A = params.A;
            elem_params.kappa = params.kappa;
            elem_params.rho = params.rho;
    end

    % 装配循环
    for elem = 1:NElem
        % 计算单元矩阵
        Ke = stiffness_func(elem_params);
        Me = mass_func(elem_params);

        % 计算节点自由度索引
        node1 = elem;
        node2 = elem + 1;

        % 节点自由度映射
        node1_dofs = (node1 - 1) * DOF_per_node + (1:DOF_per_node);
        node2_dofs = (node2 - 1) * DOF_per_node + (1:DOF_per_node);
        elem_dofs = [node1_dofs, node2_dofs];

        % 装配到全局矩阵
        K(elem_dofs, elem_dofs) = K(elem_dofs, elem_dofs) + Ke;
        M(elem_dofs, elem_dofs) = M(elem_dofs, elem_dofs) + Me;
    end
end

function Me = element(params)
% element - 轴向振动单元质量矩阵
%
% 物理模型: 一致质量矩阵 (Consistent Mass Matrix)
% 自由度: [u₁, u₂] - 轴向位移
%
% 输入:
%   params - ElementParameters对象
%
% 输出:
%   Me - 2×2 单元质量矩阵

    rho = params.rho;  % 密度
    A = params.A;      % 截面积
    L = params.L;      % 单元长度

    % 质量系数
    m_coeff = rho * A * L / 6;

    % 一致质量矩阵 (线性形状函数)
    Me = m_coeff * [2, 1; 1, 2];
end
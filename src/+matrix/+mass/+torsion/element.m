function Me = element(params)
% element - 扭转振动单元质量矩阵
%
% 物理模型: 一致质量矩阵 (Consistent Mass Matrix)
% 自由度: [φ₁, φ₂] - 扭转角
%
% 输入:
%   params - ElementParameters对象
%
% 输出:
%   Me - 2×2 单元质量矩阵

    rho = params.rho;  % 密度
    Ip = params.Ip;    % 极惯性矩
    L = params.L;      % 单元长度

    % 扭转质量系数
    m_coeff = rho * Ip * L / 6;

    % 一致质量矩阵 (线性形状函数)
    Me = m_coeff * [2, 1; 1, 2];
end
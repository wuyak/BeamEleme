function Ke = element(params)
% element - 扭转振动单元刚度矩阵
%
% 物理模型: 标准扭转杆单元 (Torsion Bar Element)
% 自由度: [φ₁, φ₂] - 扭转角
%
% 输入:
%   params - ElementParameters对象
%
% 输出:
%   Ke - 2×2 单元刚度矩阵

    G = params.G;      % 剪切模量
    J = params.J;      % 扭转惯性矩
    L = params.L;      % 单元长度

    % 扭转刚度系数
    k_torsion = G * J / L;

    % 标准扭转杆单元刚度矩阵
    Ke = k_torsion * [1, -1; -1, 1];
end
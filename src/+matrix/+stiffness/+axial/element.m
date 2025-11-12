function Ke = element(params)
% element - 轴向振动单元刚度矩阵
%
% 物理模型: 标准杆单元 (Stadium Bar Element)
% 自由度: [u₁, u₂] - 轴向位移
%
% 输入:
%   params - ElementParameters对象
%
% 输出:
%   Ke - 2×2 单元刚度矩阵

    E = params.E;      % 弹性模量
    A = params.A;      % 截面积
    L = params.L;      % 单元长度

    % 轴向刚度系数
    k_axial = E * A / L;

    % 标准杆单元刚度矩阵
    Ke = k_axial * [1, -1; -1, 1];
end
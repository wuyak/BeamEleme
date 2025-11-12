function Ke = element(params)
% element - 弯曲剪切单元刚度矩阵
%
% 物理模型: Timoshenko梁单元 (考虑剪切变形)
% 自由度: [w₁, θ₁, w₂, θ₂] - 横向位移和转角
%
% 输入:
%   params - ElementParameters对象
%
% 输出:
%   Ke - 4×4 单元刚度矩阵

    E = params.E;      % 弹性模量
    I = params.I;      % 惯性矩
    G = params.G;      % 剪切模量
    A = params.A;      % 截面积
    kappa = params.kappa;  % 剪切修正系数
    L = params.L;      % 单元长度

    % 弯曲刚度
    EI = E * I;
    % 剪切刚度
    GA = kappa * G * A;

    % Timoshenko梁单元刚度矩阵
    % 包含弯曲项和剪切项的耦合
    phi = 12 * EI / (GA * L^2);  % 剪切影响因子

    Ke = zeros(4, 4);

    % 弯曲刚度项
    Ke(1,1) =  12 * EI / (L^3 * (1 + phi));
    Ke(1,2) =  6 * EI / (L^2 * (1 + phi));
    Ke(1,3) = -12 * EI / (L^3 * (1 + phi));
    Ke(1,4) =  6 * EI / (L^2 * (1 + phi));

    Ke(2,1) =  6 * EI / (L^2 * (1 + phi));
    Ke(2,2) =  (4 + phi) * EI / (L * (1 + phi));
    Ke(2,3) = -6 * EI / (L^2 * (1 + phi));
    Ke(2,4) =  (2 - phi) * EI / (L * (1 + phi));

    Ke(3,1) = -12 * EI / (L^3 * (1 + phi));
    Ke(3,2) = -6 * EI / (L^2 * (1 + phi));
    Ke(3,3) =  12 * EI / (L^3 * (1 + phi));
    Ke(3,4) = -6 * EI / (L^2 * (1 + phi));

    Ke(4,1) =  6 * EI / (L^2 * (1 + phi));
    Ke(4,2) =  (2 - phi) * EI / (L * (1 + phi));
    Ke(4,3) = -6 * EI / (L^2 * (1 + phi));
    Ke(4,4) =  (4 + phi) * EI / (L * (1 + phi));
end
function Me = element(params)
% element - 弯曲剪切单元质量矩阵 (Friedman 1993)
%
% 物理模型: Timoshenko梁一致质量矩阵（考虑剪切变形影响）
% 自由度: [w₁, θ₁, w₂, θ₂] - 横向位移和转角
%
% 理论基础:
%   Friedman(1993) Eq.(14a, 14b) - 精确一致质量矩阵
%   Eq.(14a): 平动惯性（ρA项）- 含phi耦合
%   Eq.(14b): 转动惯性（ρI项）- 含phi耦合
%
%
% 参考文献:
%   Friedman, Z., & Kosmatka, J. B. (1993).
%   An improved two-node Timoshenko beam finite element.
%   Computers & structures, 47(3), 473-481.
%   DOI: 10.1016/0045-7949(93)90243-7
%
% 输入:
%   params - ElementParameters对象，包含:
%       .rho   - 密度 [kg/m³]
%       .A     - 截面积 [m²]
%       .I     - 惯性矩 [m⁴]
%       .E     - 弹性模量 [Pa]
%       .G     - 剪切模量 [Pa]
%       .kappa - 剪切修正系数（无量纲）
%       .L     - 单元长度 [m]
%
% 输出:
%   Me - 4×4 单元质量矩阵

    % 提取参数
    rho = params.rho;
    A = params.A;
    I = params.I;
    E = params.E;
    G = params.G;
    kappa = params.kappa;
    L = params.L;

    % 单位长度质量和转动惯量
    rhoA = rho * A;
    rhoI = rho * I;

    % 计算刚度参数
    EI = E * I;
    kGA = kappa * G * A;  % κGA（剪切刚度）

    % Friedman(1993) 无量纲参数
    phi = 12 * EI / (L^2 * kGA);

    % Eq.(14a): 平动惯性（ρA项）- 含phi耦合
    MrhoA = (rhoA * L / (210 * (1+phi)^2)) * [
        (70*phi^2 + 147*phi + 78),    (35*phi^2 + 77*phi + 44)*L/4, ...
        (35*phi^2 + 63*phi + 27),    -(35*phi^2 + 63*phi + 26)*L/4;

        (35*phi^2 + 77*phi + 44)*L/4,  (7*phi^2 + 14*phi + 8)*L^2/4, ...
        (35*phi^2 + 63*phi + 26)*L/4,  -(7*phi^2 + 14*phi + 6)*L^2/4;

        (35*phi^2 + 63*phi + 27),     (35*phi^2 + 63*phi + 26)*L/4, ...
        (70*phi^2 + 147*phi + 78),   -(35*phi^2 + 77*phi + 44)*L/4;

        -(35*phi^2 + 63*phi + 26)*L/4,  -(7*phi^2 + 14*phi + 6)*L^2/4, ...
        -(35*phi^2 + 77*phi + 44)*L/4,  (7*phi^2 + 14*phi + 8)*L^2/4
    ];

    % Eq.(14b): 转动惯性（ρI项）- 含phi耦合
    MrhoI = (rhoI / (30 * (1+phi)^2 * L)) * [
        36,                    -(15*phi - 3)*L,           -36,                   -(15*phi - 3)*L;
        -(15*phi - 3)*L,      (10*phi^2 + 5*phi + 4)*L^2, (15*phi - 3)*L,       (5*phi^2 - 5*phi - 1)*L^2;
        -36,                   (15*phi - 3)*L,            36,                    (15*phi - 3)*L;
        -(15*phi - 3)*L,      (5*phi^2 - 5*phi - 1)*L^2, (15*phi - 3)*L,       (10*phi^2 + 5*phi + 4)*L^2
    ];

    % 总质量矩阵 = 平动惯性 + 转动惯性
    Me = MrhoA + MrhoI;

end

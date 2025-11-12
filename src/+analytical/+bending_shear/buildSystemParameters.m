function sys = buildSystemParameters(params)
% buildSystemParameters - 构建Timoshenko梁的系统参数结构
%
% 输入:
%   params - BeamParameters对象
%
% 输出:
%   sys - 系统参数结构体，包含：
%     .par.L, .par.A, .par.I, .par.rho  - 物理参数
%     .par.alpha, .par.beta, .par.gamma - 无量纲参数
%     .par.T                             - 时间尺度因子
%     .numStable, .maxMode, .tol         - 算法选项
%
% 理论依据:
%   Khasawneh & Segalman (2019) 定义的无量纲参数:
%   alpha = A*L^2/I
%   beta  = ((A*G*k)/(E*I))*L^2
%   gamma = beta/alpha = G*k/E
%   T     = L*sqrt(rho/(G*k))

% 提取物理参数
E = params.E;           % 杨氏模量 (Pa)
G = params.G;           % 剪切模量 (Pa)
rho = params.rho;       % 密度 (kg/m³)
L = params.L;           % 梁长 (m)
A = params.A;           % 截面积 (m²)
I = params.I;           % 惯性矩 (m⁴)
kappa = params.kappa;     % 剪切修正系数

% 计算无量纲参数
alphaN = A * L^2 / I;
betaN = ((A * G * kappa) / (E * I)) * L^2;
gammaN = betaN / alphaN;
T = L * sqrt(rho / (G * kappa));

% 构建参数结构
sys.par.L = L;
sys.par.A = A;
sys.par.I = I;
sys.par.rho = rho;
sys.par.alpha = alphaN;
sys.par.beta = betaN;
sys.par.gamma = gammaN;
sys.par.T = T;

% 算法选项
sys.numStable = true;           % 使用数值稳定的表达式
sys.maxMode = 100000;           % 最大模态值
sys.tol = 1e-8;                 % 容差
sys.lambda_scale = 1;           % 不缩放

end

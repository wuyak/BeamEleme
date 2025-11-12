function ana_result = solve(params, BC, num_modes)
% solve - 计算扭转振动的解析解（频率+模态振型）
%
% 输入:
%   params    - BeamParameters对象
%   BC        - 边界条件 ('PP', 'CF', 'CC', 'FF')
%   num_modes - 模态数
%
% 输出:
%   ana_result - 解析解结构体
%     .freq_Hz     - 频率 [Hz] (n×1)
%     .omega       - 角频率 [rad/s] (n×1)
%     .getModes    - 函数句柄 @(x) 返回(length(x) × n)的模态矩阵（连续函数）
%     .BC          - 边界条件
%
% 示例:
%   ana_result = analytical.torsion.solve(params, 'CF', 10);
%   x = linspace(0, params.L, 100);
%   phi_modes = ana_result.getModes(x);  % 100×10矩阵

% 验证输入
assert(length(BC) == 2, '边界条件必须是2字符');
BC = upper(BC);

% 计算频率
freq_Hz = analytical.torsion.getFrequencies(params, BC, num_modes);
omega = freq_Hz * 2 * pi;

% 确保为列向量
freq_Hz = freq_Hz(:);
omega = omega(:);

% 创建模态获取函数
L = params.L;
getModes = @(x) analytical.torsion.getModes(x, L, BC, num_modes);

% 打包返回（标量结构体）
ana_result.freq_Hz = freq_Hz;
ana_result.omega = omega;
ana_result.getModes = getModes;
ana_result.BC = BC;

% 添加可视化所需的标准字段
ana_result.matrix_method = 'torsion';
ana_result.source = 'analytical';
ana_result.params = params;

end

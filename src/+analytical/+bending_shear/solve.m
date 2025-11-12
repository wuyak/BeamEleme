function ana_result = solve(params, BC, num_modes)
% solve - 计算Timoshenko梁的解析解（频率+模态振型）
%
%
% 输入:
%   params    - BeamParameters对象
%   BC        - 边界条件 ('PP', 'CF', 'CC', etc.)
%   num_modes - 模态数
%
% 输出:
%   ana_result - 解析解结构体
%     .freq_Hz     - 频率 [Hz] (n×1)
%     .omega       - 角频率 [rad/s] (n×1)
%     .lambda      - 无量纲特征值 (n×1)
%     .modes       - chebfun模态表示（连续函数）
%     .sys         - 系统参数（用于内部）
%     .BC          - 边界条件
%
% 示例:
%   mat = parameters.MaterialLibrary.Steel();
%   sec = parameters.SectionLibrary.Rectangular(0.2, 0.1);
%   params = parameters.BeamParameters(mat, sec);
%   ana_result = analytical.bending_shear.solve(params, 'PP', 10);
%
%   频率: ana_result.freq_Hz
%   模态: ana_result.modes{区域编号}(z坐标, 模态编号)

% 验证输入
assert(length(BC) == 2, '边界条件必须是2字符');
BC = upper(BC);

% 确保chebfun可用
ensureChebfunAvailable();

% 1. 构建系统参数
sys = analytical.bending_shear.buildSystemParameters(params);

% 警告：短粗梁的已知限制
if sys.par.alpha < 100
    warning('analytical:bending_shear:shortBeam', ...
        ['短粗梁(α=%.1f < 100): 解析解在 λ≈α 附近可能有缺失模态。\n' ...
         '这是特征值解析算法的缺陷，具体分析请查看 docs/analytical_solution_limitations.md'], ...
        sys.par.alpha);
end

% 2. 获取特征值的chebfun表示
eigVal_chebfun = analytical.bending_shear.getEigenvalues(sys, BC);

% 3. 求解lambda
lambda_vals = analytical.bending_shear.solveLambdaFromChebfun(eigVal_chebfun, sys.par.alpha);

% 4. 限制模态数
if length(lambda_vals) < num_modes
    warning('只找到 %d 个模态，少于请求的 %d 个', length(lambda_vals), num_modes);
    num_modes = length(lambda_vals);
end
lambda_vals = lambda_vals(1:num_modes);

% 5. 计算频率
% 根据Khasawneh论文: lambda = (omega * T)^2
omega = sqrt(lambda_vals) / sys.par.T;
freq_Hz = omega / (2 * pi);

% 确保为列向量
freq_Hz = freq_Hz(:);
omega = omega(:);

% 6. 计算模态振型
lambda_cell = analytical.bending_shear.packLambdaForModes(lambda_vals, sys.par.alpha);
cheb_modes = analytical.bending_shear.getModes(sys, BC, lambda_cell);

% 7. 打包返回（使用标量结构体）
ana_result.freq_Hz = freq_Hz;
ana_result.omega = omega;
ana_result.lambda = lambda_vals;
ana_result.modes = cheb_modes;
ana_result.sys = sys;
ana_result.BC = BC;

% 添加可视化所需的标准字段
ana_result.matrix_method = 'bending_shear';
ana_result.source = 'analytical';
ana_result.params = params;

end

%% 工具函数
function ensureChebfunAvailable()
if ~exist('chebfun', 'file')
    chebfun_path = fullfile(fileparts(mfilename('fullpath')), '..', '..', '..', 'chebfun');
    if exist(chebfun_path, 'dir')
        addpath(chebfun_path);
        chebfunroot;
        if ~exist('chebfun', 'file')
            error('无法加载 chebfun。请检查路径：%s', chebfun_path);
        end
    else
        error('需要chebfun工具箱。请安装或检查路径：%s', chebfun_path);
    end
end
end

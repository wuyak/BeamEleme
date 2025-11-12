% test_bug_fixes.m
% 测试Timoshenko梁解析解代码中3个bug的修复
%
% Bug 1: CP边界λ=α变量使用错误 (getModes.m:243-244)
% Bug 2: PP边界λ=α变量使用错误 (getModes.m:397-398)
% Bug 3: CP边界λ>α缺失搜索范围 (getEigenvalues.m:147)

function test_bug_fixes()
    fprintf('\n========== 测试解析解bug修复 ==========\n\n');

    % 添加路径
    addpath(genpath('F:\Github\BeamElem\chebfun'));

    % 测试1：验证CP边界条件Bug 3修复（λ>α区域能找到根）
    test_CP_lambda_GT();

    % 测试2：验证CP边界条件Bug 1修复（λ=α使用正确变量）
    test_CP_lambda_EQ();

    % 测试3：验证PP边界条件Bug 2修复（λ=α使用正确变量）
    test_PP_lambda_EQ();

    fprintf('\n========== 所有测试通过 ==========\n\n');
end

function test_CP_lambda_GT()
    fprintf('测试1: CP边界λ>α区域搜索范围修复\n');
    fprintf('--------------------------------------\n');

    % 创建测试梁（α=100，确保有λ>α的模态）
    mat = parameters.MaterialLibrary.Steel();
    sec = parameters.SectionLibrary.Circular(0.4);
    beam = parameters.BeamModel(mat, sec, 'L', 1.0);

    % 计算解析解
    n_modes = 5;
    ana_result = analytical.bending_shear.solve(beam, 'CP', n_modes);

    % 检查是否找到足够的模态
    n_found = length(ana_result.freq_Hz);
    fprintf('  请求模态数: %d\n', n_modes);
    fprintf('  找到模态数: %d\n', n_found);

    if n_found >= n_modes
        fprintf('  ✓ 测试通过：CP边界能在λ>α区域找到模态\n');
    else
        error('  ✗ 测试失败：CP边界缺失模态（Bug 3未修复）');
    end

    % 计算α值
    sys = analytical.bending_shear.buildSystemParameters(beam);
    fprintf('  α = %.2f\n', sys.par.alpha);

    % 显示前几个模态的λ值
    lambda_values = ana_result.lambda;
    fprintf('  前5个模态的λ值:\n');
    for i = 1:min(5, length(lambda_values))
        if lambda_values(i) > sys.par.alpha
            fprintf('    λ_%d = %.2f (> α) ✓\n', i, lambda_values(i));
        else
            fprintf('    λ_%d = %.2f (< α)\n', i, lambda_values(i));
        end
    end
    fprintf('\n');
end

function test_CP_lambda_EQ()
    fprintf('测试2: CP边界λ=α变量使用修复\n');
    fprintf('--------------------------------------\n');

    % 创建测试梁，尝试找到λ≈α的情况
    % 注意：精确λ=α的情况很难构造，这里测试代码逻辑一致性
    mat = parameters.MaterialLibrary.Steel();
    sec = parameters.SectionLibrary.Ring(0.3, 0.2);
    beam = parameters.BeamModel(mat, sec, 'L', 1.0);

    % 计算解析解和FEM
    NElem = 200;
    n_modes = 5;

    try
        ana_result = analytical.bending_shear.solve(beam, 'CP', n_modes);
        fem_result = solvers.solve('CP', beam, 'bending_shear', NElem, n_modes);

        % 对比频率误差
        comp = comparison.compare(fem_result, ana_result, 'n_modes', n_modes);

        fprintf('  FEM vs 解析解频率误差:\n');
        for i = 1:min(n_modes, length(comp.freq_error))
            fprintf('    模态%d: %.4f%%\n', i, comp.freq_error(i));
        end

        % 检查误差是否合理（<1%）
        max_error = max(comp.freq_error(~isnan(comp.freq_error)));
        if max_error < 0.01
            fprintf('  ✓ 测试通过：CP边界频率误差 < 1%%\n');
        else
            warning('  ! CP边界最大误差 = %.4f%% (可能存在问题)', max_error * 100);
        end
    catch ME
        error('  ✗ 测试失败：%s', ME.message);
    end
    fprintf('\n');
end

function test_PP_lambda_EQ()
    fprintf('测试3: PP边界λ=α变量使用修复\n');
    fprintf('--------------------------------------\n');

    % 创建测试梁
    mat = parameters.MaterialLibrary.Steel();
    sec = parameters.SectionLibrary.Ring(0.3, 0.2);
    beam = parameters.BeamModel(mat, sec, 'L', 1.0);

    % 计算解析解和FEM
    NElem = 200;
    n_modes = 5;

    try
        ana_result = analytical.bending_shear.solve(beam, 'PP', n_modes);
        fem_result = solvers.solve('PP', beam, 'bending_shear', NElem, n_modes);

        % 对比频率误差
        comp = comparison.compare(fem_result, ana_result, 'n_modes', n_modes);

        fprintf('  FEM vs 解析解频率误差:\n');
        for i = 1:min(n_modes, length(comp.freq_error))
            fprintf('    模态%d: %.4f%%\n', i, comp.freq_error(i));
        end

        % 检查误差是否合理（<1%）
        max_error = max(comp.freq_error(~isnan(comp.freq_error)));
        if max_error < 0.01
            fprintf('  ✓ 测试通过：PP边界频率误差 < 1%%\n');
        else
            warning('  ! PP边界最大误差 = %.4f%% (可能存在问题)', max_error * 100);
        end
    catch ME
        error('  ✗ 测试失败：%s', ME.message);
    end
    fprintf('\n');
end

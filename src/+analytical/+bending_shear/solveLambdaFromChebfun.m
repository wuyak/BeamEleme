function lambda_vals = solveLambdaFromChebfun(eigVal_chebfun, alphaN)
% solveLambdaFromChebfun - 从chebfun特征值表示求解lambda值
%
% 输入:
%   eigVal_chebfun - 特征值的chebfun表示（3个区域的cell数组）
%                    {1}: lambda < alpha 区域
%                    {2}: lambda = alpha 特殊点
%                    {3}: lambda > alpha 区域
%   alphaN         - 无量纲参数 alpha = A*L^2/I
%
% 输出:
%   lambda_vals - 特征值数组（已排序、去重、过滤接近0的值）
%
% 算法:
%   1. 处理三个区域的特征值
%   2. 对区域1和3使用chebfun求根
%   3. 严格过滤确保根在正确区域内
%   4. 排序、去重、去除刚体模态

lambda_vals = [];

% 处理三个区域: lambda < alpha, lambda = alpha, lambda > alpha
for i = 1:3
    if ~isempty(eigVal_chebfun{i})
        if i == 2
            % lambda = alpha 的特殊情况
            if eigVal_chebfun{i} > 0
                lambda_vals = [lambda_vals; eigVal_chebfun{i}];
            end
        else
            % lambda < alpha (i=1) 或 lambda > alpha (i=3)
            try
                roots_i = roots(eigVal_chebfun{i});
                if ~isempty(roots_i)
                    % 严格过滤：确保根在正确的区域内
                    if i == 1
                        % 区域 1: 只保留 lambda < alpha 的根
                        roots_i = roots_i(roots_i < alphaN);
                    elseif i == 3
                        % 区域 3: 只保留 lambda > alpha 的根
                        roots_i = roots_i(roots_i > alphaN);
                    end
                    lambda_vals = [lambda_vals; roots_i];
                end
            catch
                % 某些边界条件可能在某些区域没有根
                continue;
            end
        end
    end
end

% 排序
lambda_vals = sort(lambda_vals);

% 去重：使用相对容差去除重复的特征值
% 某些边界条件（如CP）会在不同区域返回相同的根
lambda_vals = uniquetol(lambda_vals, 1e-8, 'DataScale', 1);

% 去除接近0的值 (刚体模态)
lambda_vals = lambda_vals(lambda_vals > 1e-6);

end

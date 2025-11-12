function lambda_cell = packLambdaForModes(lambda_vals, alphaN)
% packLambdaForModes - 将lambda数组打包为getTimoshenkoModes需要的格式
%
% getTimoshenkoModes期望的输入格式是cell数组：
%   {1}: lambda < alpha 的值
%   {2}: lambda = alpha 的值（如果有）
%   {3}: lambda > alpha 的值
%
% 输入:
%   lambda_vals - lambda值数组（已排序）
%   alphaN      - 无量纲参数 alpha = A*L^2/I
%
% 输出:
%   lambda_cell - 分区的cell数组 {region1, region2, region3}

% 初始化三个区域
lambda_cell = cell(3, 1);

% 设置容差
tol = 1e-8;

% 分配lambda到对应区域
region1 = lambda_vals(lambda_vals < alphaN - tol);
region2 = lambda_vals(abs(lambda_vals - alphaN) <= tol);
region3 = lambda_vals(lambda_vals > alphaN + tol);

% 打包为cell数组
lambda_cell{1} = region1;
lambda_cell{2} = region2;
lambda_cell{3} = region3;

end

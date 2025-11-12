clc; clear;

% 添加src到路径
script_dir = fileparts(mfilename('fullpath'));         % src/+workflow/tests/convergence
tests_dir = fileparts(script_dir);                     % src/+workflow/tests
workflow_dir = fileparts(tests_dir);                   % src/+workflow
src_dir = fileparts(workflow_dir);                     % src
addpath(src_dir);

% 创建输出目录
output_dir = fullfile(script_dir, 'output');
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

fprintf('========================================\n');
fprintf('收敛性分析测试\n');
fprintf('========================================\n');
fprintf('结果将保存到: %s\n', output_dir);
fprintf('========================================\n\n');

%% 配置参数
% 3个截面
sections = {
    struct('type', 'rectangular', 'b', 0.2, 'h', 0.1, 'name', '矩形(0.2×0.1)', 'filename', 'rectangular_0.2x0.1'), ...
    struct('type', 'circular', 'd', 0.2, 'name', '圆形(d=0.2)', 'filename', 'circular_d0.2'), ...
    struct('type', 'ring', 'd_out', 0.3, 'd_in', 0.2, 'name', '圆环(D=0.3,d=0.2)', 'filename', 'ring_D0.3_d0.2')
};

% 完整测试条件
BCs = {'PP', 'CF', 'CP', 'PF', 'CC', 'FF', 'PR', 'CR', 'RR', 'RF'};  % 10种边界条件
NElem_list = [64, 128, 256, 512];  % 4种单元数
n_modes = 50;  % 固定模态数
L = 1.0;  % 固定长度

% 物理模态类型
physics_types = {'axial', 'torsion', 'bending_shear'};

%% 数据存储
% 结构：error_data{section_idx}{physics_idx}(BC_idx, NElem_idx)
error_data = cell(length(sections), length(physics_types));
for i = 1:length(sections)
    for j = 1:length(physics_types)
        error_data{i}{j} = nan(length(BCs), length(NElem_list));
    end
end

%% 执行收敛性测试
mat = parameters.MaterialLibrary.Steel();

for sec_idx = 1:length(sections)
    sec_config = sections{sec_idx};

    % 创建截面
    if strcmp(sec_config.type, 'rectangular')
        sec = parameters.SectionLibrary.Rectangular(sec_config.b, sec_config.h);
    elseif strcmp(sec_config.type, 'circular')
        sec = parameters.SectionLibrary.Circular(sec_config.d);
    elseif strcmp(sec_config.type, 'ring')
        sec = parameters.SectionLibrary.Ring(sec_config.d_out, sec_config.d_in);
    end

    % 创建梁模型
    beam = parameters.BeamModel(mat, sec, 'L', L);

    fprintf('\n【截面 %d/%d】%s\n', sec_idx, length(sections), sec_config.name);
    fprintf('========================================\n');

    for bc_idx = 1:length(BCs)
        BC = BCs{bc_idx};

        for ne_idx = 1:length(NElem_list)
            NElem = NElem_list(ne_idx);

            fprintf('  测试 %s, NElem=%d ... ', BC, NElem);

            try
                % 调用workflow求解（静默模式）
                wf = workflow.solve(BC, beam, NElem, n_modes, false,"convergence");

                % 提取各物理模态的最大误差
                for phy_idx = 1:length(physics_types)
                    phy = physics_types{phy_idx};
                    comp = wf.(phy).comp;

                    % 使用最大频率误差作为收敛指标
                    error_data{sec_idx}{phy_idx}(bc_idx, ne_idx) = comp.summary.max_error;
                end

                fprintf('✓\n');
            catch ME
                fprintf('✗ 失败: %s\n', ME.message);
            end
        end
    end
end

fprintf('\n========================================\n');
fprintf('数据收集完成，开始绘图...\n');
fprintf('========================================\n\n');

%% 绘制收敛曲线
physics_names = {'轴向', '扭转', '弯曲剪切'};
colors = lines(length(BCs));

for sec_idx = 1:length(sections)
    sec_config = sections{sec_idx};

    figure('Position', [100 + (sec_idx-1)*100, 100, 1200, 400]);

    for phy_idx = 1:length(physics_types)
        subplot(1, 3, phy_idx);
        hold on;
        grid on;

        % 绘制每个边界条件的收敛曲线
        for bc_idx = 1:length(BCs)
            errors = error_data{sec_idx}{phy_idx}(bc_idx, :);
            plot(NElem_list, errors, 'o-', 'LineWidth', 2, ...
                 'MarkerSize', 8, 'Color', colors(bc_idx, :), ...
                 'DisplayName', BCs{bc_idx});
        end

        set(gca, 'XScale', 'log', 'YScale', 'log');
        xlabel('单元数 NElem');
        ylabel('最大频率误差 (%)');
        title(sprintf('%s (模态数=%d)', physics_names{phy_idx}, n_modes));
        legend('Location', 'best');
    end

    sgtitle(sprintf('收敛性分析 - %s', sec_config.name), 'FontSize', 14, 'FontWeight', 'bold');

    % 保存图片到output目录
    png_path = fullfile(output_dir, sprintf('convergence_%s.png', sec_config.filename));
    saveas(gcf, png_path);
    fprintf('  ✓ 已保存图片: %s\n', png_path);
end

fprintf('\n');

%% 保存误差表
fprintf('保存误差数据...\n');

for sec_idx = 1:length(sections)
    sec_config = sections{sec_idx};

    % 保存到output目录
    filename = fullfile(output_dir, sprintf('error_table_%s.csv', sec_config.filename));

    % 打开文件
    fid = fopen(filename, 'w');

    % 写入表头
    fprintf(fid, '截面,"%s"\n', sec_config.name);
    fprintf(fid, '\n');

    for phy_idx = 1:length(physics_types)
        fprintf(fid, '%s - 最大频率误差 (%%)\n', physics_names{phy_idx});
        fprintf(fid, 'BC,');
        for ne_idx = 1:length(NElem_list)
            fprintf(fid, '%d', NElem_list(ne_idx));
            if ne_idx < length(NElem_list)
                fprintf(fid, ',');
            end
        end
        fprintf(fid, '\n');

        for bc_idx = 1:length(BCs)
            fprintf(fid, '%s,', BCs{bc_idx});
            for ne_idx = 1:length(NElem_list)
                err = error_data{sec_idx}{phy_idx}(bc_idx, ne_idx);
                fprintf(fid, '%.4f', err);
                if ne_idx < length(NElem_list)
                    fprintf(fid, ',');
                end
            end
            fprintf(fid, '\n');
        end
        fprintf(fid, '\n');
    end

    fclose(fid);
    fprintf('  ✓ 已保存CSV: %s\n', filename);
end

fprintf('\n========================================\n');
fprintf('收敛性分析完成！\n');
fprintf('========================================\n');

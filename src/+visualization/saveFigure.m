function saveFigure(fig_handle, filename, varargin)
% saveFigure - 保存图形到文件
%
% 用法:
%   visualization.saveFigure(fig, 'myplot')
%   visualization.saveFigure(fig, 'myplot', 'format', {'pdf', 'png'})
%   visualization.saveFigure(gcf, 'figure1', 'format', 'pdf', 'dpi', 300)
%
% 输入:
%   fig_handle - figure句柄
%   filename   - 文件名（不含扩展名）
%
% Name-Value参数:
%   'format'      - 格式: 'pdf'|'png'|'eps'|'svg'|'fig' 或 cell数组（默认'pdf'）
%   'dpi'         - PNG的DPI（默认300）
%   'transparent' - 是否透明背景（默认false）
%   'path'        - 保存路径（默认当前目录）
%
% 示例:
%   fig = visualization.plot(result, 'type', 'mode');
%   visualization.saveFigure(fig, 'mode_shape', 'format', {'pdf', 'png'});

% 解析参数
p = inputParser;
addRequired(p, 'fig_handle');
addRequired(p, 'filename', @ischar);
addParameter(p, 'format', 'pdf', @(x) ischar(x) || iscell(x));
addParameter(p, 'dpi', 300, @isnumeric);
addParameter(p, 'transparent', false, @islogical);
addParameter(p, 'path', '.', @ischar);
parse(p, fig_handle, filename, varargin{:});

formats = p.Results.format;
dpi_val = p.Results.dpi;
transparent = p.Results.transparent;
save_path = p.Results.path;

% 确保format是cell数组
if ischar(formats)
    formats = {formats};
end

% 创建保存路径
if ~exist(save_path, 'dir')
    mkdir(save_path);
end

% 保存每种格式
for i = 1:length(formats)
    fmt = lower(formats{i});
    output_file = fullfile(save_path, [filename '.' fmt]);

    switch fmt
        case 'pdf'
            % PDF矢量格式
            exportgraphics(fig_handle, output_file, ...
                'ContentType', 'vector', ...
                'BackgroundColor', getBackgroundColor(transparent));
            fprintf('✓ Saved: %s\n', output_file);

        case 'png'
            % PNG位图格式
            exportgraphics(fig_handle, output_file, ...
                'Resolution', dpi_val, ...
                'BackgroundColor', getBackgroundColor(transparent));
            fprintf('✓ Saved: %s (DPI=%d)\n', output_file, dpi_val);

        case 'eps'
            % EPS矢量格式
            print(fig_handle, output_file, '-depsc', '-tiff');
            fprintf('✓ Saved: %s\n', output_file);

        case 'svg'
            % SVG矢量格式
            saveas(fig_handle, output_file, 'svg');
            fprintf('✓ Saved: %s\n', output_file);

        case 'fig'
            % MATLAB figure格式
            savefig(fig_handle, output_file);
            fprintf('✓ Saved: %s\n', output_file);

        otherwise
            warning('不支持的格式: %s', fmt);
    end
end

end


function bg_color = getBackgroundColor(transparent)
% getBackgroundColor - 获取背景颜色

if transparent
    bg_color = 'none';
else
    bg_color = 'white';
end

end

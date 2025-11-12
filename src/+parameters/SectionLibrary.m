classdef SectionLibrary
% SectionLibrary - 标准截面库
%
% 提供常用截面形状，简化截面定义
%
% 用法：
%   sec = parameters.SectionLibrary.Rectangular(0.2, 0.1);
%   sec = parameters.SectionLibrary.Circular(0.15);
%   sec = parameters.SectionLibrary.Ring(0.2, 0.18);
%
% 截面结构体包含：
%   .type - 截面类型 ('rectangular', 'circular', 'ring')
%   .dims - 尺寸参数 (具体内容取决于type)
%
% 注意：
%   - 截面库只负责定义几何形状
%   - 截面属性计算(A, I, J, k)由SectionCalculator完成

    methods (Static)

        function sec = Rectangular(b, h)
            % 矩形截面
            %
            % 输入:
            %   b - 宽度 (z方向, m)，要求 b >= h
            %   h - 高度 (y方向, 弯曲方向, m)
            %
            % 约定: b >= h (Roark和Hutchinson公式要求)
            %
            % 示例:
            %   sec = parameters.SectionLibrary.Rectangular(0.2, 0.1);

            assert(b > 0 && h > 0, '截面尺寸必须大于0');
            assert(b >= h, '矩形截面要求 b >= h（宽度>=高度），当前 b=%.3f, h=%.3f', b, h);

            sec = struct(...
                'type', 'rectangular', ...
                'dims', struct('b', b, 'h', h));
        end


        function sec = Circular(d)
            % 圆形截面
            %
            % 输入:
            %   d - 直径 (m)
            %
            % 示例:
            %   sec = parameters.SectionLibrary.Circular(0.15);

            assert(d > 0, '直径必须大于0');

            sec = struct(...
                'type', 'circular', ...
                'dims', struct('d', d));
        end


        function sec = Ring(d_out, d_in)
            % 环形截面（空心圆）
            %
            % 输入:
            %   d_out - 外径 (m)
            %   d_in  - 内径 (m)
            %
            % 示例:
            %   sec = parameters.SectionLibrary.Ring(0.2, 0.18);

            assert(d_out > 0 && d_in > 0, '截面尺寸必须大于0');
            assert(d_out > d_in, '外径必须大于内径，当前 d_out=%.3f, d_in=%.3f', d_out, d_in);

            sec = struct(...
                'type', 'ring', ...
                'dims', struct('d_out', d_out, 'd_in', d_in));
        end


        function sec = IBeam(designation)
            % 工字钢截面（预留，待实现）
            %
            % 输入:
            %   designation - 型号，如 'I20a'
            %
            % TODO: 添加标准工字钢尺寸表
            %
            % 示例:
            %   sec = parameters.SectionLibrary.IBeam('I20a');

            error('IBeam暂未实现。请使用 Rectangular, Circular 或 Ring');
        end


        function sec = HBeam(designation)
            % H型钢截面（预留，待实现）
            %
            % 输入:
            %   designation - 型号，如 'H200x100x5x8'
            %
            % TODO: 添加标准H型钢尺寸表
            %
            % 示例:
            %   sec = parameters.SectionLibrary.HBeam('H200x100x5x8');

            error('HBeam暂未实现。请使用 Rectangular, Circular 或 Ring');
        end


        function listSections()
            % 列出截面库使用示例
            %
            % 示例:
            %   parameters.SectionLibrary.listSections();

            fprintf('\n========================================\n');
            fprintf('截面库\n');
            fprintf('========================================\n\n');

            fprintf('【矩形截面】\n');
            fprintf('  sec = parameters.SectionLibrary.Rectangular(b, h);\n');
            fprintf('  示例: Rectangular(0.2, 0.1)  // 200mm×100mm\n');
            fprintf('  约定: b >= h (宽度 >= 高度)\n\n');

            fprintf('【圆形截面】\n');
            fprintf('  sec = parameters.SectionLibrary.Circular(d);\n');
            fprintf('  示例: Circular(0.15)  // 直径150mm\n\n');

            fprintf('【环形截面】\n');
            fprintf('  sec = parameters.SectionLibrary.Ring(d_out, d_in);\n');
            fprintf('  示例: Ring(0.2, 0.18)  // 外径200mm, 内径180mm\n\n');

            fprintf('【预留接口】\n');
            fprintf('  IBeam(''I20a'')    - 工字钢 (待实现)\n');
            fprintf('  HBeam(''H200...'')  - H型钢 (待实现)\n\n');

            fprintf('========================================\n');
            fprintf('用法示例:\n');
            fprintf('  mat = parameters.MaterialLibrary.Steel();\n');
            fprintf('  sec = parameters.SectionLibrary.Rectangular(0.2, 0.1);\n');
            fprintf('  params = parameters.BeamParameters(mat, sec, ''NElem'', 100);\n');
            fprintf('========================================\n\n');
        end


        function printSection(sec)
            % 打印截面信息
            %
            % 输入:
            %   sec - 截面结构体
            %
            % 示例:
            %   sec = parameters.SectionLibrary.Rectangular(0.2, 0.1);
            %   parameters.SectionLibrary.printSection(sec);

            fprintf('截面类型: %s\n', sec.type);

            switch sec.type
                case 'rectangular'
                    fprintf('  宽度 b = %.3f m\n', sec.dims.b);
                    fprintf('  高度 h = %.3f m\n', sec.dims.h);

                case 'circular'
                    fprintf('  直径 d = %.3f m\n', sec.dims.d);

                case 'ring'
                    fprintf('  外径 d_out = %.3f m\n', sec.dims.d_out);
                    fprintf('  内径 d_in  = %.3f m\n', sec.dims.d_in);
                    fprintf('  壁厚 t = %.3f m\n', (sec.dims.d_out - sec.dims.d_in)/2);

                otherwise
                    fprintf('  未知截面类型\n');
            end
        end

    end
end

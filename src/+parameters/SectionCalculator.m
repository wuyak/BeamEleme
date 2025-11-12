classdef SectionCalculator
% SectionCalculator - 截面属性计算器
%
% 功能：
%   从截面几何(section)和材料属性(material)计算截面属性(A, I, J, Ip, k)
%
% 设计原则：
%   - 单一职责：只负责计算，不负责存储
%   - 数据驱动：用映射表替代if/else
%   - 可测试：纯函数，无副作用
%
% 用法：
%   mat = parameters.MaterialLibrary.Steel();
%   sec = parameters.SectionLibrary.Rectangular(0.2, 0.1);
%   props = parameters.SectionCalculator.calculate(sec, mat);
%
% 输出结构体:
%   props.A   - 截面积 (m²)
%   props.I   - 惯性矩 Iz (m⁴)
%   props.J   - 抗扭惯性矩 (m⁴)
%   props.Ip  - 极面积惯性矩 (m⁴)
%   props.kappa - 剪切修正系数 (无量纲)
%   props.G   - 剪切模量 (Pa)
%
% 参考文献：
%   - Roark's Formulas for Stress and Strain (截面参数)
%   - Hutchinson 2001 (剪切修正系数)

    methods (Static)

        function props = calculate(section, material)
            % 计算截面属性（主入口）
            %
            % 输入:
            %   section  - 截面结构体 (来自SectionLibrary)
            %   material - 材料结构体 (来自MaterialLibrary)
            %
            % 输出:
            %   props - 截面属性结构体

            % 计算剪切模量
            G = material.E / (2 * (1 + material.nu));

            % 根据截面类型调用相应计算函数
            switch section.type
                case 'rectangular'
                    props = parameters.SectionCalculator.rectangular(section.dims, material.nu);

                case 'circular'
                    props = parameters.SectionCalculator.circular(section.dims, material.nu);

                case 'ring'
                    props = parameters.SectionCalculator.ring(section.dims, material.nu);

                otherwise
                    error('SectionCalculator:UnsupportedType', ...
                          '未支持的截面类型: %s', section.type);
            end

            % 添加剪切模量
            props.G = G;
        end


        function props = rectangular(dims, nu)
            % 矩形截面属性计算
            %
            % 输入:
            %   dims - struct('b', ..., 'h', ...)
            %   nu   - 泊松比
            %
            % 参考: Roark's Ch.10 Table 7, Hutchinson 2001 式(46)(47)

            b = dims.b;
            h = dims.h;

            % 截面积和惯性矩 (Roark's Table A.1, Case 2)
            A = b * h;
            I = b * h^3 / 12;  % Iz = b*h³/12，绕z轴，用于xy平面弯曲

            % 极面积惯性矩 (Ip = Iy + Iz)
            Iy = h * b^3 / 12;  % 绕y轴
            Iz = b * h^3 / 12;  % 绕z轴
            Ip = Iy + Iz;

            % 抗扭惯性矩 (Roark's Table 10.7, Case 4)
            % J = ab³[16/3 - 3.36(b/a)(1 - b⁴/12a⁴)]  for a ≥ b
            % 这里 a=b（长边）, b_short=h（短边）
            J = b * h^3 * (16/3 - 3.36*(h/b)*(1 - h^4/(12*b^4)));

            % 剪切修正系数 (Hutchinson 2001, 式46+47)
            a = b / 2;
            b_half = h / 2;

            % 计算C4（式46）- 向量化无穷级数
            C4_part1 = (4/45) * a^3 * b_half * (-12*a^2 - 15*nu*a^2 + 5*nu*b_half^2);

            n = 1:50;
            npi = n * pi;
            terms = (16*nu^2*b_half^5) ./ ((npi).^5 * (1+nu)) .* ...
                    (npi*a/b_half) .* tanh(npi*a/b_half);
            C4_part2 = nu * sum(terms);
            C4 = C4_part1 + C4_part2;

            % 计算k（式47）
            kappa = -2*(1+nu) / ((9/(4*a^5*b_half))*C4 + nu*(1 - b_half^2/a^2));

            % 封装返回
            props = struct('A', A, 'I', I, 'J', J, 'Ip', Ip, 'kappa', kappa);
        end


        function props = circular(dims, nu)
            % 圆形截面属性计算
            %
            % 输入:
            %   dims - struct('d', ...)
            %   nu   - 泊松比
            %
            % 参考: Roark's Table A.1 Case 15, Hutchinson 2001 式(43)

            d = dims.d;

            % 截面积和惯性矩 (Roark's Table A.1, Case 15)
            A = pi * d^2 / 4;
            I = pi * d^4 / 64;

            % 极面积惯性矩 (Ip = 2*I for circular section)
            Ip = pi * d^4 / 32;

            % 抗扭惯性矩 (Roark's Table 10.7, Case 1)
            J = pi * d^4 / 32;

            % 剪切修正系数 (Hutchinson 2001, 式43)
            kappa = 6*(1+nu)^2 / (7 + 12*nu + 4*nu^2);

            % 封装返回
            props = struct('A', A, 'I', I, 'J', J, 'Ip', Ip, 'kappa', kappa);
        end


        function props = ring(dims, nu)
            % 环形截面属性计算（圆形空心截面）
            %
            % 输入:
            %   dims - struct('d_out', ..., 'd_in', ...)
            %   nu   - 泊松比
            %
            % 参考: Roark's Table A.1 Case 16, Hutchinson 2001 式(44)(53)

            d_out = dims.d_out;
            d_in = dims.d_in;

            % 截面积和惯性矩 (Roark's Table A.1, Case 16)
            A = pi * (d_out^2 - d_in^2) / 4;
            I = pi * (d_out^4 - d_in^4) / 64;

            % 极面积惯性矩 (Ip = 2*I for circular section)
            Ip = pi * (d_out^4 - d_in^4) / 32;

            % 抗扭惯性矩 (Roark's Table 10.7, Case 10)
            J = pi * (d_out^4 - d_in^4) / 32;

            % 剪切修正系数 (Hutchinson 2001)
            alpha = d_in / d_out;

            if alpha < 0.001
                % 实心圆（α → 0）：退化为圆形截面公式（式43）
                kappa = 6*(1+nu)^2 / (7 + 12*nu + 4*nu^2);
            elseif alpha > 0.95
                % 薄壁管（α → 1）：Hutchinson式(53)
                kappa = (1 + nu) / (2 + nu);
            else
                % 一般环形截面：Hutchinson式(44)完整公式
                numerator = 6 * (alpha^2 + 1)^2 * (1 + nu)^2;
                c1 = 7 + 12*nu + 4*nu^2;
                c2 = 34 + 48*nu + 16*nu^2;
                denominator = c1*alpha^4 + c2*alpha^2 + c1;
                kappa = numerator / denominator;
            end

            % 封装返回
            props = struct('A', A, 'I', I, 'J', J, 'Ip', Ip, 'kappa', kappa);
        end


        function printProps(props)
            % 打印截面属性
            %
            % 输入:
            %   props - 截面属性结构体
            %
            % 示例:
            %   parameters.SectionCalculator.printProps(props);

            fprintf('截面属性:\n');
            fprintf('  A  = %.4e m²\n', props.A);
            fprintf('  I  = %.4e m⁴\n', props.I);
            fprintf('  J  = %.4e m⁴\n', props.J);
            fprintf('  Ip = %.4e m⁴\n', props.Ip);
            fprintf('  κ  = %.4f\n', props.kappa);
            fprintf('  G  = %.2e Pa\n', props.G);
        end

    end
end

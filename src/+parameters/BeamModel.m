classdef BeamModel < handle
% BeamModel - 梁物理模型类
%
% 用法（简洁！）:
%   mat = parameters.MaterialLibrary.Steel();
%   sec = parameters.SectionLibrary.Rectangular(0.2, 0.1);
%   beam = parameters.BeamModel(mat, sec, 'L', 1.0);
%
% 核心设计：
%   - 只包含梁的物理属性（材料、截面、几何）
%   - 不包含离散化参数（NElem）或求解配置（num_modes）
%   - 材料和截面从库中选取（MaterialLibrary + SectionLibrary）
%   - 截面计算委托给SectionCalculator（分离职责）

    properties
        % === 材料属性 ===
        E           % 弹性模量 (Pa)
        G           % 剪切模量 (Pa)
        rho         % 密度 (kg/m³)
        nu          % 泊松比

        % === 截面属性 ===
        section_type % 'rectangular', 'circular', 'ring'
        A           % 截面积 (m²)
        I           % 惯性矩 Iz (m⁴)
        J           % 抗扭惯性矩 (m⁴)
        Ip          % 极面积惯性矩 (m⁴)
        kappa       % 剪切修正系数

        % 截面尺寸（仅用于记录，不参与计算）
        dims        % 尺寸结构体 (与SectionLibrary一致)

        % === 几何参数 ===
        L           % 梁长 (m)

        % === 刚度参数 ===
        EI          % 弯曲刚度
        EA          % 轴向刚度
        GJ          % 扭转刚度

        % === 惯性参数 ===
        rhoA        % 单位长度质量
        rhoI        % 单位长度转动惯量

        % === 波速参数 ===
        c           % 纵波波速
        c_bending   % 弯曲波速
    end

    methods
        function obj = BeamModel(material, section, varargin)
            % 构造函数
            %
            % 输入:
            %   material - 材料结构体 (来自MaterialLibrary)
            %   section  - 截面结构体 (来自SectionLibrary)
            %   varargin - Name-Value参数对
            %              'L', 1.0          % 梁长
            %
            % 示例:
            %   mat = parameters.MaterialLibrary.Steel();
            %   sec = parameters.SectionLibrary.Rectangular(0.2, 0.1);
            %   beam = parameters.BeamModel(mat, sec, 'L', 1.0);

            % 设置材料属性
            obj.E = material.E;
            obj.rho = material.rho;
            obj.nu = material.nu;

            % 计算截面属性（委托给SectionCalculator）
            props = parameters.SectionCalculator.calculate(section, material);
            obj.A = props.A;
            obj.I = props.I;
            obj.J = props.J;
            obj.Ip = props.Ip;
            obj.kappa = props.kappa;
            obj.G = props.G;

            % 记录截面类型和尺寸
            obj.section_type = section.type;
            obj.dims = section.dims;

            % 设置默认参数
            obj.L = 1.0;

            % 应用用户自定义参数
            obj.setOptions(varargin{:});

            % 更新派生参数
            obj.updateDerivedParameters();
        end


        function setOptions(obj, varargin)
            % 设置Name-Value参数
            %
            % 示例:
            %   beam.setOptions('L', 2.0);

            for i = 1:2:length(varargin)
                prop_name = varargin{i};
                prop_value = varargin{i+1};

                if isprop(obj, prop_name)
                    obj.(prop_name) = prop_value;
                else
                    warning('BeamModel:UnknownProperty', ...
                            '未知属性: %s', prop_name);
                end
            end

            % 更新派生参数
            obj.updateDerivedParameters();
        end


        function printParameters(obj)
            % 打印参数摘要

            fprintf('\n========================================\n');
            fprintf('梁参数摘要\n');
            fprintf('========================================\n\n');

            fprintf('【材料】\n');
            fprintf('  E   = %.2e Pa\n', obj.E);
            fprintf('  G   = %.2e Pa\n', obj.G);
            fprintf('  ρ   = %.0f kg/m³\n', obj.rho);
            fprintf('  ν   = %.2f\n\n', obj.nu);

            fprintf('【截面】 %s\n', obj.section_type);
            switch obj.section_type
                case 'rectangular'
                    fprintf('  b   = %.3f m (宽度)\n', obj.dims.b);
                    fprintf('  h   = %.3f m (高度)\n', obj.dims.h);
                case 'circular'
                    fprintf('  d   = %.3f m (直径)\n', obj.dims.d);
                case 'ring'
                    fprintf('  d_out = %.3f m (外径)\n', obj.dims.d_out);
                    fprintf('  d_in  = %.3f m (内径)\n', obj.dims.d_in);
            end
            fprintf('  A   = %.4e m²\n', obj.A);
            fprintf('  I   = %.4e m⁴\n', obj.I);
            fprintf('  J   = %.4e m⁴\n', obj.J);
            fprintf('  κ   = %.4f\n\n', obj.kappa);

            fprintf('【几何】\n');
            fprintf('  L   = %.2f m\n\n', obj.L);

            fprintf('\n========================================\n\n');
        end


        function validateParameters(obj)
            % 验证参数有效性

            assert(obj.E > 0, '弹性模量必须大于0');
            assert(obj.G > 0, '剪切模量必须大于0');
            assert(obj.rho > 0, '密度必须大于0');
            assert(obj.nu >= 0 && obj.nu < 0.5, '泊松比必须在 [0, 0.5) 范围');
            assert(obj.L > 0, '梁长必须大于0');
            assert(obj.A > 0, '截面积必须大于0');
            assert(obj.I > 0, '惯性矩必须大于0');
            assert(obj.J > 0, '抗扭惯性矩必须大于0');
            assert(obj.Ip > 0, '极面积惯性矩必须大于0');
            assert(obj.kappa > 0, '剪切修正系数必须大于0');
        end


        function hash = toHash(obj)
            % toHash - 生成参数对象的唯一哈希值
            %
            % 只包含物理参数（梁的固有属性）：
            %   - 材料: E, G, rho
            %   - 截面: A, I, J, kappa
            %   - 几何: L
            %
            % 不包含离散化参数（这些会显式放在缓存键中）：
            %   - NElem, num_modes (在缓存键中显式表示为 N{NElem}_M{num_modes})
            %
            % 输出:
            %   hash - 8位16进制字符串 (例如: 'A3F2B91C')
            %
            % 用途:
            %   生成缓存键的一部分，表示梁的物理属性
            %   同一根梁，不同离散化参数 → 相同哈希值

            % 打包物理参数为字符串（不含 NElem, num_modes）
            key_params = sprintf('E%.6e_G%.6e_rho%.6e_A%.6e_I%.6e_J%.6e_L%.6f_kappa%.6f', ...
                obj.E, obj.G, obj.rho, obj.A, obj.I, obj.J, obj.L, obj.kappa);

            % 计算哈希（使用 MATLAB 内置 Java String.hashCode）
            hash_value = java.lang.String(key_params).hashCode();

            % 转为8位16进制字符串（处理负数）
            if hash_value < 0
                hash_value = hash_value + 2^32;
            end
            hash = upper(dec2hex(hash_value, 8));
        end


        function updateDerivedParameters(obj)
            % 更新派生参数
            %
            % 包括: 刚度、惯性、波速
            % 注: node_type 硬编码为 2（两节点线性单元）

            % 刚度参数
            obj.EI = obj.E * obj.I;
            obj.EA = obj.E * obj.A;
            obj.GJ = obj.G * obj.J;

            % 惯性参数
            obj.rhoA = obj.rho * obj.A;
            obj.rhoI = obj.rho * obj.I;

            % 波速参数
            obj.c = sqrt(obj.E / obj.rho);
            obj.c_bending = sqrt(obj.EI / obj.rhoA);
        end
    end
end

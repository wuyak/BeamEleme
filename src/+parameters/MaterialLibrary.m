classdef MaterialLibrary
% MaterialLibrary - 标准材料库
%
% 提供常用工程材料的参数，避免手动输入易出错
%
% 用法：
%   mat = parameters.MaterialLibrary.Steel();
%   mat = parameters.MaterialLibrary.Aluminum();
%   mat = parameters.MaterialLibrary.Custom('Titanium', 110e9, 4500, 0.34);
%
% 材料结构体包含：
%   .name - 材料名称
%   .E    - 弹性模量 (Pa)
%   .rho  - 密度 (kg/m³)
%   .nu   - 泊松比 (无量纲)
%
% 支持添加自定义材料：
%   直接在此文件中添加新的静态方法

    methods (Static)

        function mat = Steel(grade)
            % 结构钢（默认Q235）
            %
            % 输入:
            %   grade - (可选) 钢材牌号
            %           'Q235' (默认) - 普通结构钢
            %           'Q345' - 低合金高强度钢
            %           'Q420' - 高强度结构钢
            %
            % 示例:
            %   mat = parameters.MaterialLibrary.Steel();       % Q235
            %   mat = parameters.MaterialLibrary.Steel('Q345');

            if nargin < 1
                grade = 'Q235';
            end

            switch grade
                case 'Q235'
                    mat = struct(...
                        'name', 'Q235结构钢', ...
                        'E', 200e9, ...      % Pa
                        'rho', 7850, ...     % kg/m³
                        'nu', 0.3);

                case 'Q345'
                    mat = struct(...
                        'name', 'Q345低合金钢', ...
                        'E', 206e9, ...
                        'rho', 7850, ...
                        'nu', 0.3);

                case 'Q420'
                    mat = struct(...
                        'name', 'Q420高强钢', ...
                        'E', 206e9, ...
                        'rho', 7850, ...
                        'nu', 0.29);

                otherwise
                    error('未知钢材牌号: %s (支持: Q235, Q345, Q420)', grade);
            end
        end


        function mat = Aluminum(alloy)
            % 铝合金
            %
            % 输入:
            %   alloy - (可选) 合金牌号
            %           '6061' (默认) - 结构用铝合金
            %           '7075' - 高强度铝合金
            %
            % 示例:
            %   mat = parameters.MaterialLibrary.Aluminum();        % 6061
            %   mat = parameters.MaterialLibrary.Aluminum('7075');

            if nargin < 1
                alloy = '6061';
            end

            switch alloy
                case '6061'
                    mat = struct(...
                        'name', '6061铝合金', ...
                        'E', 69e9, ...       % Pa
                        'rho', 2700, ...     % kg/m³
                        'nu', 0.33);

                case '7075'
                    mat = struct(...
                        'name', '7075铝合金', ...
                        'E', 72e9, ...
                        'rho', 2810, ...
                        'nu', 0.33);

                otherwise
                    error('未知铝合金牌号: %s (支持: 6061, 7075)', alloy);
            end
        end


        function mat = Concrete(grade)
            % 混凝土
            %
            % 输入:
            %   grade - (可选) 强度等级
            %           'C30' (默认)
            %           'C40', 'C50'
            %
            % 示例:
            %   mat = parameters.MaterialLibrary.Concrete();      % C30
            %   mat = parameters.MaterialLibrary.Concrete('C40');

            if nargin < 1
                grade = 'C30';
            end

            switch grade
                case 'C30'
                    mat = struct(...
                        'name', 'C30混凝土', ...
                        'E', 30e9, ...       % Pa
                        'rho', 2400, ...     % kg/m³
                        'nu', 0.2);

                case 'C40'
                    mat = struct(...
                        'name', 'C40混凝土', ...
                        'E', 32.5e9, ...
                        'rho', 2400, ...
                        'nu', 0.2);

                case 'C50'
                    mat = struct(...
                        'name', 'C50混凝土', ...
                        'E', 34.5e9, ...
                        'rho', 2400, ...
                        'nu', 0.2);

                otherwise
                    error('未知混凝土等级: %s (支持: C30, C40, C50)', grade);
            end
        end


        function mat = Timber(species)
            % 木材
            %
            % 输入:
            %   species - (可选) 木材种类
            %             'pine' (默认) - 松木
            %             'oak' - 橡木
            %
            % 示例:
            %   mat = parameters.MaterialLibrary.Timber();        % 松木
            %   mat = parameters.MaterialLibrary.Timber('oak');

            if nargin < 1
                species = 'pine';
            end

            switch species
                case 'pine'
                    mat = struct(...
                        'name', '松木', ...
                        'E', 10e9, ...       % Pa (沿纹理)
                        'rho', 500, ...      % kg/m³
                        'nu', 0.37);

                case 'oak'
                    mat = struct(...
                        'name', '橡木', ...
                        'E', 12e9, ...
                        'rho', 720, ...
                        'nu', 0.37);

                otherwise
                    error('未知木材种类: %s (支持: pine, oak)', species);
            end
        end


        function mat = Custom(name, E, rho, nu)
            % 自定义材料
            %
            % 输入:
            %   name - 材料名称
            %   E    - 弹性模量 (Pa)
            %   rho  - 密度 (kg/m³)
            %   nu   - 泊松比 (无量纲)
            %
            % 示例:
            %   mat = parameters.MaterialLibrary.Custom('Titanium', 110e9, 4500, 0.34);

            assert(E > 0, '弹性模量必须大于0');
            assert(rho > 0, '密度必须大于0');
            assert(nu >= 0 && nu < 0.5, '泊松比必须在 [0, 0.5) 范围内');

            mat = struct(...
                'name', name, ...
                'E', E, ...
                'rho', rho, ...
                'nu', nu);
        end


        function listMaterials()
            % 列出所有预定义材料
            %
            % 示例:
            %   parameters.MaterialLibrary.listMaterials();

            fprintf('\n========================================\n');
            fprintf('材料库\n');
            fprintf('========================================\n\n');

            % 钢材
            fprintf('【结构钢】\n');
            materials = {'Q235', 'Q345', 'Q420'};
            for i = 1:length(materials)
                mat = parameters.MaterialLibrary.Steel(materials{i});
                fprintf('  %-8s: E=%.0e Pa, ρ=%.0f kg/m³, ν=%.2f\n', ...
                    mat.name, mat.E, mat.rho, mat.nu);
            end
            fprintf('\n');

            % 铝合金
            fprintf('【铝合金】\n');
            alloys = {'6061', '7075'};
            for i = 1:length(alloys)
                mat = parameters.MaterialLibrary.Aluminum(alloys{i});
                fprintf('  %-12s: E=%.0e Pa, ρ=%.0f kg/m³, ν=%.2f\n', ...
                    mat.name, mat.E, mat.rho, mat.nu);
            end
            fprintf('\n');

            % 混凝土
            fprintf('【混凝土】\n');
            grades = {'C30', 'C40', 'C50'};
            for i = 1:length(grades)
                mat = parameters.MaterialLibrary.Concrete(grades{i});
                fprintf('  %-8s: E=%.0e Pa, ρ=%.0f kg/m³, ν=%.2f\n', ...
                    mat.name, mat.E, mat.rho, mat.nu);
            end
            fprintf('\n');

            % 木材
            fprintf('【木材】\n');
            species = {'pine', 'oak'};
            for i = 1:length(species)
                mat = parameters.MaterialLibrary.Timber(species{i});
                fprintf('  %-4s: E=%.0e Pa, ρ=%.0f kg/m³, ν=%.2f\n', ...
                    mat.name, mat.E, mat.rho, mat.nu);
            end

            fprintf('\n========================================\n');
            fprintf('用法示例:\n');
            fprintf('  mat = parameters.MaterialLibrary.Steel();           %% Q235\n');
            fprintf('  mat = parameters.MaterialLibrary.Aluminum(''7075''); %% 高强铝合金\n');
            fprintf('  mat = parameters.MaterialLibrary.Custom(''Ti6Al4V'', 114e9, 4430, 0.34);\n');
            fprintf('========================================\n\n');
        end

    end
end


function n_rigid = getRigidModeCount(matrix_method, BC)
% getRigidModeCount - 根据物理性质和边界条件返回刚体模态数
%
% 输入:
%   matrix_method - 物理性质 ('axial', 'torsion', 'bending_shear')
%   BC            - 边界条件 (如 'FF', 'CC', 'RF', 'RR' 等)
%
% 输出:
%   n_rigid - 刚体模态数量
%
% 说明:
%   按物理性质分类：
%   - Axial + FF       → 1个刚体模态（整体平移）
%   - Torsion + FF/PF/PP     → 1个刚体模态（整体转动）
%   - BendingShear + FF → 2个刚体模态（平移+转动）
%   - BendingShear + RF/RR/PF → 1个刚体模态
%   - 其他情况         → 0个刚体模态

switch matrix_method
    case 'axial'
        % 轴向：映射后为FF边界有1个刚体模态（整体平移）
        BC_mapped = analytical.axial.mapBoundaryConditionToAxial(BC);
        if strcmp(BC_mapped, 'FF')
            n_rigid = 1;
        else
            n_rigid = 0;
        end

    case 'torsion'
        % 扭转：映射后为FF边界有1个刚体模态（整体转动）
        BC_mapped = analytical.torsion.mapBoundaryConditionToTorsion(BC);
        if strcmp(BC_mapped, 'FF')
            n_rigid = 1;
        else
            n_rigid = 0;
        end

    case 'bending_shear'
        % 弯曲剪切：FF有2个刚体模态，RF/RR有1个刚体模态
        if strcmp(BC, 'FF')
            n_rigid = 2;
        elseif ismember(BC, {'RF', 'RR', 'PF'})
            n_rigid = 1;
        else
            n_rigid = 0;
        end

    otherwise
        error('不支持的 matrix_method: %s', matrix_method);
end

end

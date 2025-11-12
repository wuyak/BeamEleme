function DOF_per_node = getDOFPerNode(matrix_method)
% getDOFPerNode - 从物理性质推断每节点自由度数
%
% 输入:
%   matrix_method - 物理类型 ('axial' | 'torsion' | 'bending_shear')
%
% 输出:
%   DOF_per_node - 每节点自由度数
%                  axial/torsion: 1
%                  bending_shear: 2
%
% 说明:
%   这是 1-to-1 映射，DOF_per_node 完全由 matrix_method 决定

switch matrix_method
    case {'axial', 'torsion'}
        DOF_per_node = 1;
    case 'bending_shear'
        DOF_per_node = 2;
    otherwise
        error('getDOFPerNode:InvalidMethod', ...
            '未知的 matrix_method: %s。仅支持 axial, torsion, bending_shear', matrix_method);
end

end

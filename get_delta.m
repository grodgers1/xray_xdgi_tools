function [ delta ] = get_delta( material, lambda )
%GET_DELTA Gives decrement of the real part of the index of refraction
%   material: e.g. 'Si' lambda: wavelength [m] 
%   NOTE: when adding a new material (say the Nth material), add the
%   material table as 'material_table_N.mat'
%   NOTE: tables should have lambda in column 1 and delta in column 2
load('materials_list.mat', 'materials_list')
table = find(materials_list == string(material));
load(['material_table_' num2str(table, '%02d') '.mat'], delta_table);
index = find(min(delta_table(:,1)-lambda));
delta = delta_table(index,2);

end


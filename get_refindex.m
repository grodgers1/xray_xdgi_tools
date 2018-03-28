function [ delta, beta ] = get_refindex( material, lambda )
%GET_REFINDEX Gives the index of refraction (1-delta+ibeta)
%   material: e.g. 'Si' 
%   lambda: wavelength [m] 
%   NOTE: when adding a new material (e.g. 'Si'), add the
%   material table as 'material_table_Si.mat'
%   NOTE: tables should have lambda in column 1, delta in column 2, beta in column 3
table = struct2cell(load(['material_table_' material '.mat']));
table = table{1};
index = find(min(table(:,1)-lambda));
delta = table(index,2);
beta = table(index,3);
end


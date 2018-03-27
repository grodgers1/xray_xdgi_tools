function [ delta, beta ] = get_refindex( material, lambda )
%GET_REFINDEX Gives the index of refraction (1-delta+ibeta)
%   material: e.g. 'Si' 
%   lambda: wavelength [m] 
%   NOTE: when adding a new material (e.g. 'Si'), add the
%   material table as 'material_table_Si.mat'
%   NOTE: tables should have lambda in column 1, delta in column 2, beta in column 3
load(['material_table_' material '.mat']);
index = find(min(au(:,1)-lambda));
delta = au(index,2);
beta = au(index,3);
end


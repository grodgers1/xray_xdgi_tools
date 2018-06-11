function [ delta, beta ] = get_refindex( material, E )
%GET_REFINDEX Gives the index of refraction (1-delta+ibeta)
%   material: e.g. 'Si' 
%   lambda: wavelength [m] 
%   NOTE: when adding a new material (e.g. 'Si'), add the
%   material table as 'material_table_Si.mat'
%   NOTE: tables should have lambda in column 1, delta in column 2, beta in column 3
material_table = struct2cell(load([material '_table.mat']));
betas = material_table{1};
deltas = material_table{2};
energies = material_table{3}*1e3;
[~,index] = min(abs(energies-E));
delta = deltas(index);
beta = betas(index);
end


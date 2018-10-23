function [design_E] = find_design_E(m1,t1)
% % Default setting is to find energy at which the m1 and t1 give pi shift

% load table
material_table = struct2cell(load([m1 '_table.mat']));
betas = material_table{1};
deltas = material_table{2};
energies = material_table{3}*1e3;
% convert energy to wavenumber
wavenumbers = 2*pi./lambda_from_E(energies);
% get phase shift (projection approximation)
shifts = t1*wavenumbers.*deltas;
% find index of energy that pi shifts
[~,index] = min(abs(shifts-pi));
% output that energy
design_E = energies(index);
end


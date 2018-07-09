%% Big simulation
% loops over the following:
%material_range = 1:3;
E_range = (25:2.5:50)*1e3;
p1_range = (3:0.25:9)*1e-6;
t1_range = (2.5:1:20)*1e-6;
d_sg1_range = (5:1:55)*1e-2;
z_range = (5:1:55)*1e-2;
%dutycycle_range = 0.3:0.1:0.7;

%%

sensitivity_mat = zeros(length(E_range),length(p1_range),length(t1_range),length(d_sg1_range),length(z_range));

n_sims = numel(sensitivity_mat);









%% Big simulation
addpath('./tables')
addpath('./functions')
clear all
% clc
%% Loops over the following:
%material_range = 1:3;
E_range = (25:5:35)*1e3;
p1_range = (3:0.5:9)*1e-6;
t1_range = (2:2:20)*1e-6;
d_sg1_range = (15:2:55)*1e-2;
z_range = (5:2:45)*1e-2;
%dutycycle_range = 0.3:0.1:0.7;

%% Initialize
sensitivity_mat = zeros(length(E_range),length(t1_range),length(p1_range),length(d_sg1_range),length(z_range));
vis_mat = zeros(length(E_range),length(t1_range),length(p1_range),length(d_sg1_range),length(z_range));
amp_mat = zeros(length(E_range),length(t1_range),length(p1_range),length(d_sg1_range),length(z_range));
counts_mat = zeros(length(E_range),length(t1_range),length(p1_range),length(d_sg1_range),length(z_range));

n_sims = numel(sensitivity_mat);

%% hardcoded values
sig_E = 4000;
n = 15;
% % get counts from amplitude and energy
% %(see ./tables/counts_measurements_nanotom.m)
% counts_energy_focus0 = @(E_mean) (0.9467*E_mean.^2 - 3.955*E_mean - 125.8);
counts_energy_focus1 = @(E_mean) (0.03238*E_mean.^2 + 22.45*E_mean - 406.5);
counts_at_60 = counts_energy_focus1(E_range/1000)/2; % divide by two for the g2

%% Run simulations
count = 1;
for e = 1:length(E_range)
    E_0 = E_range(e);
    % build spectrum
    [E_spectrum,E_x] = EspectrumGauss(E_0, sig_E, n,E_0-3*sig_E,E_0+3*sig_E);
	for t = 1:length(t1_range)
        t1 = t1_range(t);
        for p = 1:length(p1_range)
            p1 = p1_range(p);
            fprintf(['Working on sims ' num2str(count) '-' num2str(count + length(d_sg1_range)) ' of ' num2str(n_sims/length(z_range)) ' (' num2str(100*count/(n_sims/length(z_range))) ' percent) \n'])
            tic
            for dsg1 = 1:length(d_sg1_range)
                d_sg1 = d_sg1_range(dsg1);
                % build initial wavefront
                [vis,amp,~] = carpet_sim(z_range,d_sg1,p1,E_x,E_spectrum,'Au',t1);
                vis_mat(e,t,p,dsg1,:) = vis;
                amp_mat(e,t,p,dsg1,:) = amp;
                a = abs(d_sg1+z_range-0.6);
                index = find(a == min(a),1);
                counts_mat(e,t,p,dsg1,:) = amp.*(counts_at_60(e)./amp(index));
                count = count+1;
            end
            toc
        end
    end
end


%% Sensitivity calculation
% % 1/alpha_min = (pi*z/p2)*vis*sqrt(N) 
M_mat = (d_sg1_range' + z_range)./d_sg1_range';
tmp(1,:,:) = p1_range;
p2_mat(1,1,:,:,:) = permute(M_mat.*tmp,[3 1 2]);
z_mat(1,1,1,1,:) = z_range;
%tot_dist(1,1,1,:,:) = d_sg1_range' + z_range;

[ZZ,DD] = meshgrid(z_range,d_sg1_range);
totdist_mat = ZZ+DD;


sensitivity = (pi*z_mat./p2_mat).*vis_mat.*sqrt(counts_mat);
alpha_min = 1./sensitivity; % units should be [radians]

%% Visualize results
% find indices of current setup
p1_7_index = find(p1_range == 7e-06);
t1_6_index = find(t1_range == 6e-06);

vis_maps = permute(squeeze(vis_mat(:,t1_6_index,p1_7_index,:,:)),[2 3 1]);
s_maps = permute(squeeze(sensitivity(:,t1_6_index,p1_7_index,:,:)),[2 3 1]);
amin_maps = permute(squeeze(alpha_min(:,t1_6_index,p1_7_index,:,:)),[2 3 1]);

e_val = 2;

figure, imagesc(z_range, d_sg1_range, vis_maps(:,:,e_val), [0 1]), colorbar
xlabel('propagation distance [m]','interpreter','latex')
ylabel('source - $g_1$ distance [m]','interpreter','latex')
title(['Visibility for E = ' num2str(E_range(e_val)/1000) ' keV'],'interpreter','latex')
hold on, plot(z_range,0.6-z_range, 'r-')

figure, imagesc(z_range, d_sg1_range, squeeze(p2_mat(1,1,p1_7_index,:,:)), [7e-6 2.5e-5]), colorbar
xlabel('propagation distance [m]','interpreter','latex')
ylabel('source - $g_1$ distance [m]','interpreter','latex')
title(['$p_2$ for E = ' num2str(E_range(e_val)/1000) ' keV'],'interpreter','latex')

figure, imagesc(z_range, d_sg1_range, sqrt(squeeze(counts_mat(e_val,t1_6_index,p1_7_index,:,:)))), colorbar
xlabel('propagation distance [m]','interpreter','latex')
ylabel('source - $g_1$ distance [m]','interpreter','latex')
title(['$\sqrt{N}$ for E = ' num2str(E_range(e_val)/1000) ' keV'],'interpreter','latex')

figure, imagesc(z_range, d_sg1_range, s_maps(:,:,e_val)), colorbar
xlabel('propagation distance [m]','interpreter','latex')
ylabel('source - $g_1$ distance [m]','interpreter','latex')
title(['Sensitivity for E = ' num2str(E_range(e_val)/1000) ' keV'],'interpreter','latex')
hold on, plot(z_range,0.6-z_range, 'r-')

figure, imagesc(z_range, d_sg1_range, totdist_mat, [0.2 0.65]), colorbar
xlabel('propagation distance [m]','interpreter','latex')
ylabel('source - $g_1$ distance [m]','interpreter','latex')
title(['Source-detector distance [m]'],'interpreter','latex')
hold on, plot(z_range,0.6-z_range, 'r-')
    
    
    
    
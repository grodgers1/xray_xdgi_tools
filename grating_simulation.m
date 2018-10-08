%% Big simulation
addpath('./tables')
addpath('./functions')
clear all
clc
%% Loops over the following:
m1_range = 1;
materials = {'Au'};
m1 = materials{1};
E_range = (25:5:45)*1e3;
p1_range = (2:0.5:7)*1e-6;
t1_range = (2:2:12)*1e-6;
l_range = (20:1:55)*1e-2;
d_range = (5:1:40)*1e-2;
%dc_range = 0.3:0.2:0.7;
dc_range = 0.5;
%% Initialize
sensitivity_mat = zeros(length(E_range),length(t1_range),length(p1_range),length(dc_range),length(l_range),length(d_range));
vis_mat = zeros(length(E_range),length(t1_range),length(p1_range),length(dc_range),length(l_range),length(d_range));
amp_mat = zeros(length(E_range),length(t1_range),length(p1_range),length(dc_range),length(l_range),length(d_range));
counts_mat = zeros(length(E_range),length(t1_range),length(p1_range),length(dc_range),length(l_range),length(d_range));

dimension_names = {'E','t1','p1','dc','l','d'};
n_sims = numel(sensitivity_mat);
%% hardcoded values
sig_E = 4000;
n = 11;
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
            fprintf(['Working on sims ' num2str(count) '-' num2str(count + length(p1_range)*length(dc_range)) ' of ' num2str(n_sims/length(d_range)) ' (' num2str(100*count/(n_sims/length(d_range))) ' percent) \n'])
            tic
            for d = 1:length(dc_range)
            dc = dc_range(d);
                for dsg1 = 1:length(l_range)
                    d_sg1 = l_range(dsg1);
                    % build initial wavefront
                    [vis,amp,~] = carpet_sim(d_range,d_sg1,p1,dc,E_x,E_spectrum,m1,t1);
                    vis_mat(e,t,p,d,dsg1,:) = vis;
                    amp_mat(e,t,p,d,dsg1,:) = amp;
                    a = abs(d_sg1+d_range-0.6);
                    index = find(a == min(a),1);
                    counts_mat(e,t,p,d,dsg1,:) = amp.*(counts_at_60(e)./amp(index));
                    count = count+1;
                end
            end
            tmp = toc;
            fprintf(['Elapsed time ' num2str(tmp) 's. \n'])
            minrm = (n_sims/length(d_range)-count)*(tmp/(length(l_range)*length(dc_range)))/60;
            fprintf(['Estimated time remaining: ' num2str(round(minrm,1)) ' minutes. \n'])
        end    
    end
end

%% Find indices corresponding to current setup
curr_E = 30000;
curr_t1 = 6e-6;
curr_p1 = 7e-6;
curr_dc = 0.5;
curr_l = 0.296;
curr_d = 0.296;

icurr_E = find(abs(E_range-curr_E) == min(abs(E_range-curr_E)),1);
icurr_t1 = find(abs(t1_range-curr_t1) == min(abs(t1_range-curr_t1)),1);
icurr_p1 = find(abs(p1_range-curr_p1) == min(abs(p1_range-curr_p1)),1);
icurr_dc = find(abs(dc_range-curr_dc) == min(abs(dc_range-curr_dc)),1);
icurr_l = find(abs(l_range-curr_l) == min(abs(l_range-curr_l)),1);
icurr_d = find(abs(d_range-curr_d) == min(abs(d_range-curr_d)),1);
%% Sensitivity calculation
[EE,TT,PP,DCDC,LL,DD] = ndgrid(E_range,t1_range,p1_range,dc_range,l_range,d_range);
MM = (LL+DD)./LL;
p2_mat = MM.*PP/2;
totdist_mat = LL+DD;

nandist_mat = ones(size(totdist_mat));
nandist_mat(totdist_mat>0.60) = NaN; % put NaNs where total distance exceeds 60 cm

snr_alpha = (DD.*vis_mat.*sqrt(counts_mat))./(p2_mat.*EE.^2); % proportional (dropped constant)
curr_snr_alpha = snr_alpha(icurr_E,icurr_t1,icurr_p1,icurr_dc,icurr_l,icurr_d);
rel_snr_alpha = snr_alpha./curr_snr_alpha;
relnan_snr_alpha = snr_alpha.*nandist_mat/curr_snr_alpha;
%% find maximum
[max_snr, imax_snr] = max(relnan_snr_alpha(:),[],'omitnan');
[imax_E, imax_t1, imax_p1, imax_dc, imax_l, imax_d] = ind2sub(size(relnan_snr_alpha),imax_snr);

fprintf(['Optimal parameters: \n'])
fprintf(['E = ' num2str(E_range(imax_E)) '\n'])
fprintf(['t1 = ' num2str(t1_range(imax_t1)) '\n'])
fprintf(['p1 = ' num2str(p1_range(imax_p1)) '\n'])
fprintf(['dc = ' num2str(dc_range(imax_dc)) '\n'])
fprintf(['l = ' num2str(l_range(imax_l)) '\n'])
fprintf(['d = ' num2str(d_range(imax_d)) '\n'])
fprintf(['SNR improvement:' num2str(rel_snr_alpha(imax_E,imax_t1,imax_p1,imax_dc,imax_l,imax_d)) '\n'])

% m1_relnan_snr_alpha = relnan_snr_alpha(1,:,:,:,:,:,:);
% [maxm1_snr, imaxm1_snr] = max(m1_relnan_snr_alpha(:),[],'omitnan');
% [~,imaxm1_E, imaxm1_t1, imaxm1_p1, imaxm1_dc, imaxm1_l, imaxm1_d] = ind2sub(size(m1_relnan_snr_alpha),imaxm1_snr);
% 
% m2_relnan_snr_alpha = relnan_snr_alpha(2,:,:,:,:,:,:);
% [maxm2_snr, imaxm2_snr] = max(m2_relnan_snr_alpha(:),[],'omitnan');
% [~,imaxm2_E, imaxm2_t1, imaxm2_p1, imaxm2_dc, imaxm2_l, imaxm2_d] = ind2sub(size(m2_relnan_snr_alpha),imaxm2_snr);
%% Visualize results
figure, imagesc(d_range, l_range, squeeze(vis_mat(icurr_E,icurr_t1,icurr_p1,icurr_dc,:,:)), [0 1]), colorbar
xlabel('propagation distance [m]','interpreter','latex')
ylabel('source - $g_1$ distance [m]','interpreter','latex')
title(['Visibility for (E = ' num2str(E_range(icurr_E)) ', $t_1$ = ' ...
    num2str(t1_range(icurr_t1)) ', $p_1$ = ' num2str(p1_range(icurr_p1)) ...
    ', $dc$ = ' num2str(dc_range(icurr_dc)) ')'],'interpreter','latex')
hold on, plot(d_range,curr_l+curr_d-d_range, 'r-')
plot(curr_d,curr_l,'ro','MarkerSize',20)

figure, imagesc(d_range, l_range, squeeze(rel_snr_alpha(icurr_E,icurr_t1,icurr_p1,icurr_dc,:,:)),[0 max_snr]), colorbar
xlabel('propagation distance [m]','interpreter','latex')
ylabel('source - $g_1$ distance [m]','interpreter','latex')
title(['Relative $SNR_{\alpha}$  for (E = ' num2str(E_range(icurr_E)) ', $t_1$ = ' ...
    num2str(t1_range(icurr_t1)) ', $p_1$ = ' num2str(p1_range(icurr_p1)) ...
    ', $dc$ = ' num2str(dc_range(icurr_dc)) ')'],'interpreter','latex')
hold on, plot(d_range,0.6-d_range, 'r-')
plot(d_range(icurr_d),l_range(icurr_l),'ro','MarkerSize',20)

figure, imagesc(d_range, l_range, squeeze(rel_snr_alpha(imax_E,imax_t1,imax_p1,imax_dc,:,:)),[0 max_snr]), colorbar
xlabel('propagation distance [m]','interpreter','latex')
ylabel('source - $g_1$ distance [m]','interpreter','latex')
title(['Maxmum $SNR_{\alpha}$, (E = ' num2str(E_range(imax_E)) ...
    ', $t_1$ = ' num2str(t1_range(imax_t1)) ', $p_1$ = ' num2str(p1_range(imax_p1))...
    ', $dc$ = ' num2str(dc_range(imax_dc)) ')'],'interpreter','latex')
hold on, plot(d_range,0.6-d_range, 'r-')
plot(d_range(imax_d),l_range(imax_l),'r.','MarkerSize',20)
%% less important plots

%% find snr improvement of a given setup
try_E = 25000;
try_t1 = 5e-6;
try_p1 = 4e-6;
try_dc = 0.5;
try_l = 0.50;
try_d = 0.10;

itry_E = find(abs(E_range-try_E) == min(abs(E_range-try_E)),1);
itry_t1 = find(abs(t1_range-try_t1) == min(abs(t1_range-try_t1)),1);
itry_p1 = find(abs(p1_range-try_p1) == min(abs(p1_range-try_p1)),1);
itry_dc = find(abs(dc_range-try_dc) == min(abs(dc_range-try_dc)),1);
itry_l = find(abs(l_range-try_l) == min(abs(l_range-try_l)),1);
itry_d = find(abs(d_range-try_d) == min(abs(d_range-try_d)),1);

try_snr = rel_snr_alpha(itry_E,itry_t1,itry_p1,itry_dc,itry_l,itry_d);

fprintf(['Test setup: \n'])
fprintf(['E = ' num2str(E_range(itry_E)) '\n'])
fprintf(['t1 = ' num2str(t1_range(itry_t1)) '\n'])
fprintf(['p1 = ' num2str(p1_range(itry_p1)) '\n'])
fprintf(['dc = ' num2str(dc_range(itry_dc)) '\n'])
fprintf(['l = ' num2str(l_range(itry_l)) '\n'])
fprintf(['d = ' num2str(d_range(itry_d)) '\n'])
fprintf(['SNR improvement:' num2str(try_snr) '\n'])
%% vismaps
figure, imagesc(t1_range,p1_range,squeeze(vis_mat(itry_E,:,:,1,itry_l,itry_d)), [0 1])
xlabel('thickness'), ylabel('period')
title('Visibility')
colorbar

figure, imagesc(t1_range,p1_range,squeeze(rel_snr_alpha(itry_E,:,:,1,itry_l,itry_d)), [0 max_snr])
xlabel('thickness'), ylabel('period')
title('SNR')
colorbar
%% Get carpets for setups
sim_d_range = (5:0.1:30)*1e-2;
inewcurr_d = find(abs(sim_d_range-curr_d) == min(abs(sim_d_range-curr_d)),1);
inewmax_d = find(abs(sim_d_range-d_range(imax_d)) == min(abs(sim_d_range-d_range(imax_d))),1);


[E_spectrum,E_x] = EspectrumGauss(E_range(icurr_E), sig_E, n,E_range(icurr_E)-3*sig_E,E_range(icurr_E)+3*sig_E);
[vis_curr,amp_curr,carpet_curr] = carpet_sim(sim_d_range,l_range(icurr_l),p1_range(icurr_p1),dc_range(icurr_dc),E_x,E_spectrum,m1,t1_range(icurr_t1));
                    
[E_spectrum,E_x] = EspectrumGauss(E_range(imax_E), sig_E, n,E_range(imax_E)-3*sig_E,E_range(imax_E)+3*sig_E);
[vis_max,amp_max,carpet_max] = carpet_sim(sim_d_range,l_range(imax_l),p1_range(imax_p1),dc_range(imax_dc),E_x,E_spectrum,m1,t1_range(imax_t1));

y_axis_curr = (-(size(carpet_curr,1)/2):(size(carpet_curr,1)/2)-1)*1e6/(50/p1_range(icurr_p1));
figure, imagesc(sim_d_range,y_axis_curr,carpet_curr,[0 2.5]), colormap gray
hold on
vline(sim_d_range(inewcurr_d),'r-')
ylim([-12.5 12.5])
xlabel('$d$ [m]','Interpreter','latex')
ylabel('[um]','Interpreter','latex')
title('Current setup','Interpreter','latex')

y_axis_max = (-(size(carpet_max,1)/2):(size(carpet_max,1)/2)-1)*1e6/(50/p1_range(imax_p1));
figure, imagesc(sim_d_range,y_axis_max,carpet_max,[0 2.5]), colormap gray
hold on
vline(sim_d_range(inewmax_d),'r-')
ylim([-12.5 12.5])
xlabel('$d$ [m]','Interpreter','latex')
ylabel('[um]','Interpreter','latex')
title('Optimal setup','Interpreter','latex')

prof_curr = carpet_curr(:,inewcurr_d);
prof_max = carpet_max(:,inewmax_d);
figure, plot(y_axis_curr,prof_curr)
hold on, title('Current')
ylim([0 1.4])
figure, plot(y_axis_max,prof_max)
hold on, title('Optimized')
ylim([0 1.4])




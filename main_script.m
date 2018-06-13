
clear all
addpath('./tables/')
addpath('./functions/')


%% Settings
% % Define first grating
p1_range = (1:0.25:10)*1e-06; % [m] period of the above mentioned pattern
z = 0.0025:0.0025:0.4; % propagation distances [m]
steps = 20; % phase steps
detector_pixsize = 100e-6;
periods = 3;
m1 = 'Au'; % First material (e.g. 'Au')
m2 = 'Au'; % Second material (e.g. 'Au')
t1 = 6e-06; % thickness of first material [m]
t2 = 0; % thickness of second material [m]

% % Define source
%source_size = 2e-06; % FWHM of source [m]
source_size = 0;
%d_sg1 = 0.294; % source to g1 distance [m]
d_sg1 = 1.5;
% % Define spectrum
%source_spectrum = 'gaussian';
source_spectrum = 'monochromatic';
E_0 = 30000;
sig_E = 4000;
n = 25;
E_min = 15000;
E_max = 45000;
if strcmp(source_spectrum, 'gaussian') == 1 % Gaussian source:
[E_spectrum,E_x] = EspectrumGauss(E_0, sig_E, n,E_min,E_max);
else % Monochromatic source
E_x = E_0;
E_spectrum = 1;
end
figure, plot(E_x, E_spectrum, '.-'), xlabel('Energy [eV]')
ylabel('Intensity'), title('Energy Spectrum')

% % Simulation options
%wf_method = 'magnification'
wf_method = 'spherical wavefront';
x_pixels = 750; % pixels per period for simulation
reptimes = 15; % repeat g1 so that large z aren't smeared

N = x_pixels*reptimes; 
lambda = lambda_from_E(E_x)';
k = 2*pi./lambda';

vis_map = zeros(length(z),length(p1_range));
amp_map = zeros(length(z),length(p1_range));
vis_map1 = zeros(length(z),length(p1_range));
amp_map1 = zeros(length(z),length(p1_range));
curves_mat = zeros(steps,length(z),length(p1_range));
carpets_mat = zeros(N,length(z),length(p1_range));
for p = 1:length(p1_range)   
    tic
    p1 = p1_range(p);

    pixsize = p1/x_pixels;
    x = (-(N/2):(N/2-1))*pixsize; % real space coordinates
    
    %% Create initial wavefront and grating
    g1_pattern = zeros(length(x),length(E_x));
    wf1  = zeros(length(x),length(E_x));
    for e = 1:length(E_x)
        [delta1,beta1] = get_refindex(m1, E_x(e));
        [delta2,beta2] = get_refindex(m2, E_x(e));
        tmp1 = E_spectrum(e)*exp(-1i*(delta1-1i*beta1)*k(e)*t1);
        tmp2 = E_spectrum(e)*exp(-1i*(delta2-1i*beta2)*k(e)*t2);
        temp = [tmp1*ones(1,round(x_pixels/2)) tmp2*ones(1,round(x_pixels/2))]';
        temp = repmat(temp,reptimes,1);
        g1_pattern(:,e) = temp;
        if strcmp(wf_method,'spherical wavefront') == 1
            wf1(:,e)  = temp.*exp( 1i.*k(e).*sqrt(d_sg1.^2+x.^2) )'./sqrt(d_sg1.^2+x.^2)';
        else
            wf1(:,e)  = temp;
        end
    end

    %% Propagate wavefront
    I_full = wf_propagate(wf1,z,E_x,pixsize,source_size,d_sg1, 'method','spherical wavefront');
    %% phase stepping
    vis = zeros(1,length(z)); amp = zeros(1,length(z));
    curves = zeros(steps,length(z));
    for dis = 1:length(z)
        M = (d_sg1+z(dis))/d_sg1;
        p2 = p1*M/2;
        [stepping_curve]=phase_stepping(I_full(:,dis),pixsize,detector_pixsize,p2,0.5,steps,periods);
        ft_sc = fft(stepping_curve);
        vis(dis) = 2*abs(ft_sc(1+periods))./abs(ft_sc(1));
        amp(dis) = mean(stepping_curve);
        curves(:,dis) = stepping_curve;
    end
    vis_map(:,p) = vis; amp_map(:,p) = amp;
    curves_mat(:,:,p) = curves;
    carpets_mat(:,:,p) = I_full;
    toc
end

% [PP,ZZ] = meshgrid(p1_range,z); 
% relative_sensitivity = (ZZ.*vis_map.*sqrt(amp_map))./PP;

%% Calculating Talbot distances
D_1 = (p1_range.^2)/(8*lambda_from_E(E_0));
d_1 = D_1.*d_sg1./(d_sg1-D_1);
d_1(end) = d_1(end-1);
index = find(p1_range == 7e-06);
d_1(index)

%% Plotting
xtc = round(size(vis_map,1)/10);
ytc = round(size(vis_map,2)/10);

figure, imagesc(carpets_mat(:,:,index)), colormap gray
xticks(xtc:xtc:length(z))
xticklabels(z(xtc:xtc:length(z)))
yticks([])
xlabel('inter-grating distance [m]')
ylabel('position (p1 = 7 um)')

figure, imagesc(vis_map')
colorbar
yticks(ytc:ytc:size(vis_map,2))
yticklabels(p1_range(ytc:ytc:size(vis_map,2))*1e6)
xticks(xtc:xtc:size(vis_map,1))
xticklabels(z(xtc:xtc:size(vis_map,1)))
xlabel('inter-grating distance [m]')
ylabel('p1 [um]')
title('visibility')
hold on
plot(fliplr(D_1)*length(z)/(z(end)-z(1)),fliplr(1:length(p1_range)),'r-')
plot(fliplr(3*D_1)*length(z)/(z(end)-z(1)),fliplr(1:length(p1_range)),'r-')
plot(fliplr(d_1)*length(z)/(z(end)-z(1)),fliplr(1:length(p1_range)),'r-')

figure, imagesc(amp_map')
colorbar
yticks(ytc:ytc:size(amp_map,2))
yticklabels(p1_range(ytc:ytc:size(amp_map,2))*1e6)
xticks(xtc:xtc:size(amp_map,1))
xticklabels(z(xtc:xtc:size(amp_map,1)))
xlabel('inter-grating distance [m]')
ylabel('p1 [um]')
title('amplitude')

p2_mat = p1_range'.*((d_sg1+z)/d_sg1)/2; % M*p1/2 [m]
figure, imagesc(p2_mat)
colorbar
yticks(ytc:ytc:size(p2_mat,1))
yticklabels(p1_range(ytc:ytc:size(p2_mat,1))*1e6)
xticks(xtc:xtc:size(p2_mat,2))
xticklabels(z(xtc:xtc:size(p2_mat,2)))
xlabel('inter-grating distance [m]')
ylabel('p1 [um]')
title('p2 [um]')





















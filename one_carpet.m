%% Simulate one Talbot carpet
clear all
clc
addpath('./tables/')
addpath('./functions/')
addpath('./images/')

troubleshooting = 0;
%% Settings (frequently changed)
% % Propagation distances
z = 0.005:0.005:0.30; % propagation distances [m]
% % Define first grating
p1 = 7.2e-06; % [m] period of the above mentioned pattern
m1 = 'Au'; % First material (e.g. 'Au')
t1 = 12e-06; % thickness of first material [m]
dc = 0.5; % duty cycle
% % Define source
source_size = 1.5e-06; % sigma of source [m]
% % Geometry
d_sg1 = 0.40; % source to g1 distance [m]
% % Source spectrum
E_0 = 60000;
sig_E = 20000;
%% Settings (usually unchanged)
source_spectrum = 'gaussian'; %'monochromatic';
m2 = 'Au'; % Second material (e.g. 'Au')
t2 = 0; % thickness of second material [m]
n = 40;
% % Simulation options
wf_method = 'spherical wavefront';
x_pixels = 1000; % pixels per period for simulation
reptimes = 10; % repeat g1 so that large z aren't smeared
padsize = round(x_pixels*reptimes/3);
% % options for calculating visibility
vis_method = 'direct';
steps = 15; % phase steps
detector_pixsize = 100e-6;
periods = 2;
%% Timing
tic
%% Make source spectrum
if strcmp(source_spectrum, 'gaussian') == 1 % Gaussian source:
[E_spectrum,E_x] = EspectrumGauss(E_0,sig_E,n,E_0-3*sig_E,E_0+3*sig_E);
else % Monochromatic source
E_x = E_0;
E_spectrum = 1;
end
%% Define some quantities
N = x_pixels*reptimes; 
lambda = lambda_from_E(E_x)';
k = 2*pi./lambda';
pixsize = p1/x_pixels;
x = (-(N/2):(N/2-1))*pixsize; % real space coordinates   
%% Create initial wavefront and grating    
% % % NOTE: wf1 = wf0.*g1_pattern, so it is equivalent to add the intensity
% spectrum to either wf0 or g1_pattern
wf0 = exp( 1i.*k'.*sqrt(d_sg1.^2+x.^2) )'./sqrt(d_sg1.^2+x.^2)'; % before g1
g1_pattern = build_g1(m1,m2,t1,t2,dc,E_x,E_spectrum,x_pixels,reptimes);
wf1 = wf0.*g1_pattern; % after g1, using projection approximation
wf1 = padarray(wf1,[padsize 0]); % pad for simulation to avoid edge reflection
%% Troubleshooting/visualizing initial wavefront and grating
design_E = find_design_E(m1,t1);
phase_shift_g1 = build_g1_phaseshift(m1,m2,t1,t2,dc,E_x,x_pixels,reptimes);
temp_g1_pattern = build_g1(m1,m2,t1,t2,dc,E_x,ones(size(E_spectrum)),x_pixels,reptimes);

p_0 = (k'.*sqrt(d_sg1.^2+x.^2))';
I_0 = sum(abs(wf0).^2,2);
I_1 = sum(abs(wf1(padsize+1:(end-padsize),:)).^2,2);

if troubleshooting == 1
figure, imagesc(E_x*1e-3,x*1e6,abs(wf0).^2), hold on, xlabel('Energy [keV]')
ylabel('Position [um]'), colorbar, title('Intensity of wavefront before g1')
vline(design_E*1e-3, 'r', 'Design Energy'), vline(E_0*1e-3, 'g', 'Peak energy')

% figure, imagesc(E_x*1e-3,x*1e6,p_0), hold on, xlabel('Energy [keV]')
% ylabel('Position [um]'), colorbar, title('Phase of wavefront before g1')
% vline(design_E*1e-3, 'r', 'Design Energy'), vline(E_0*1e-3, 'g', 'Peak energy')

figure, imagesc(E_x*1e-3,x*1e6,abs(wf1(padsize+1:(end-padsize),:)).^2), hold on, xlabel('Energy [keV]')
ylabel('Position [um]'), colorbar, title('Intensity of wavefront after g1')
vline(design_E*1e-3, 'r', 'Design Energy'), vline(E_0*1e-3, 'g', 'Peak energy')

figure, imagesc(E_x*1e-3,x*1e6,abs(temp_g1_pattern).^2), hold on, xlabel('Energy [keV]')
ylabel('Position [um]'), colorbar, title('Transmission of g1')
vline(design_E*1e-3, 'r', 'Design Energy'), vline(E_0*1e-3, 'g', 'Peak energy')

figure, imagesc(E_x*1e-3,x*1e6,phase_shift_g1), hold on, xlabel('Energy [keV]')
ylabel('Position [um]'), colorbar, title('Phase shift of g1')
vline(design_E*1e-3, 'r', 'Design Energy'), vline(E_0*1e-3, 'g', 'Peak energy')

figure, plot(E_x*1e-3,E_spectrum), hold on, xlabel('Energy [keV]')
ylabel('Relative intensity (0-1)'), title('Energy Spectrum')

figure, plot(x*1e6, p_0(:,find(abs(E_x-design_E) == min(abs(E_x-design_E))))-min(p_0(:,find(abs(E_x-design_E) == min(abs(E_x-design_E))))))
hold on, xlabel('Position [um]'), ylabel('Phase shift [radians]')
title('Wavefront phase at design energy before g1 due to spherical source')
end
%% Propagate wavefront
I_full = wf_propagate(wf1,z,E_x,pixsize,source_size,d_sg1, 'method','spherical wavefront');
I_full = I_full(padsize+1:(end-padsize),:);
%% Troubleshooting carpet and propagation
%% Visibility calculation
M = (d_sg1+z)/d_sg1;
curves = zeros(steps,length(z)); % initialize
if strcmp(vis_method,'phase stepping') == 1
    for dis = 1:length(z)
        p2 = p1*M(dis)/2;
        [stepping_curve]=phase_stepping(I_full(:,dis),pixsize,detector_pixsize,p2,0.5,steps,periods);
        ft_sc = fft(stepping_curve);
        vis(dis) = 2*abs(ft_sc(1+periods))./abs(ft_sc(1));
        amp(dis) = mean(stepping_curve);
        curves(:,dis) = stepping_curve;
    end
else
    for dis = 1:length(z)
        psize = M(dis)*x_pixels;
        curve = I_full(round(N/2-psize/2):round(N/2+psize/2),dis);
        ft_curve = fft(curve);
        vis(dis) = 2*abs(ft_curve(1+periods))./abs(ft_curve(1));
        amp(dis) = mean(curve);
        curves(:,dis) = interp1(1:length(curve),curve,1:length(curve)/steps:length(curve))';
    end
end
%% timing
toc
%% visualize results
y_carpet = (0:x_pixels*reptimes-1)*pixsize;

carpet = I_full((round(end/2)-2*x_pixels):(round(end/2)+2*x_pixels),:);
%y_carpet = y_carpet((round(end/2)-x_pixels):(round(end/2)+x_pixels));
y_carpet = y_carpet(1:4*x_pixels+1);
y_carpet = y_carpet-y_carpet(round(end/2));
y_carpet = y_carpet*1e6;

figure, imagesc(z,y_carpet,carpet)
hold on, colorbar, vline(0.6-d_sg1,'r','detector'), colormap gray
xlabel('d [m]'), ylabel('Position [um]')

detector_index = find(abs(z-(0.6-d_sg1)) == min(abs(z-(0.6-d_sg1))));

figure, plot(y_carpet, carpet(:,detector_index))
hold on, title('Interference pattern at detector (60 cm)')
xlabel('Position [um]'), ylabel('Intensity')
vline(0,'r'), vline(M(detector_index)*p1*1e6/2, 'r', [num2str(M(detector_index)*p1*1e6/2) ' um'])

figure, plot(z,vis)
hold on, xlabel('d [m]'), ylabel('visibility')
title(['l = ' num2str(d_sg1) ' m, p_1 = ' num2str(p1*1e6) ' um , E = ' num2str(E_0*1e-3) ' keV, \sigma_E = ' num2str(sig_E*1e-3) ' keV, s = ' num2str(source_size*1e6) ' um'])
vline(z(detector_index), 'r', '60 cm')
















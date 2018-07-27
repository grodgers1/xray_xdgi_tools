%% Simulate one Talbot carpet
clear all
clc
addpath('./tables/')
addpath('./functions/')
addpath('./images/')
%% Settings
% % Propagation distances
z = 0.00025:0.00025:0.12; % propagation distances [m]

% % Define first grating
p1 = 4e-06; % [m] period of the above mentioned pattern
m1 = 'Au'; % First material (e.g. 'Au')
m2 = 'Au'; % Second material (e.g. 'Au')
t1 = 4e-06; % thickness of first material [m]
t2 = 0; % thickness of second material [m]
dc = 0.7; % duty cycle

% % Define source
source_size = 2e-06; % FWHM of source [m]
%source_size = 0;
d_sg1 = 0.48; % source to g1 distance [m]

% % Define spectrum
source_spectrum = 'gaussian';
%source_spectrum = 'monochromatic';
E_0 = 25000;
sig_E = 4000;
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
%% Make source
if strcmp(source_spectrum, 'gaussian') == 1 % Gaussian source:
[E_spectrum,E_x] = EspectrumGauss(E_0, sig_E, n,E_0-3*sig_E,E_0+3*sig_E);
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
wf0 = exp( 1i.*k'.*sqrt(d_sg1.^2+x.^2) )'./sqrt(d_sg1.^2+x.^2)'; % before g1
g1_pattern = build_g1(m1,m2,t1,t2,dc,E_x,E_spectrum,x_pixels,reptimes);
wf1 = wf0.*g1_pattern; % after g1, using projection approximation
wf1 = padarray(wf1,[padsize 0]); % pad for simulation to avoid edge reflection
%% Propagate wavefront
I_full = wf_propagate(wf1,z,E_x,pixsize,source_size,d_sg1, 'method','spherical wavefront');
I_full = I_full(padsize:(end-padsize-1),:);
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

carpet = I_full((round(end/2)-x_pixels):(round(end/2)+x_pixels),:);
y_carpet = y_carpet((round(end/2)-x_pixels):(round(end/2)+x_pixels));



figure, imagesc(z,y_carpet,carpet, [0 1.25]), colormap gray
%figure, imagesc(z,y_carpet,carpet_monocohere, [0 25]), colormap gray

%nanotom_carpet = uint8(255*(carpet/1.25));
%monocohere_carpet = uint8((255)*(carpet_monocohere/25));

% 
% savedir = './images/';
% imwrite(nanotom_carpet,[savedir 'nanotom_carpet.tif'])
% 
% 



























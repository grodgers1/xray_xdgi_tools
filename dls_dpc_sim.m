
clear all
addpath('./tables/')
addpath('./functions/')


%% Settings
% % Define first grating
p1 = 7e-06; % [m] period of the above mentioned pattern
z = 0.0025:0.0025:0.8; % propagation distances [m]
steps = 20; % phase steps
detector_pixsize = 4.8e-6;
periods = 3;
m1 = 'Au'; % First material (e.g. 'Au')
m2 = 'Au'; % Second material (e.g. 'Au')
t1 = 70e-06; % thickness of first material [m]
t2 = 0; % thickness of second material [m]

theta = 0;

% % Define source
source_size = 300e-06; % FWHM of source [m]
d_sg1 = 250; % source to g1 distance [m]

% % Define spectrum
source_spectrum = 'monochromatic';
E_0 = 20000;

% % Simulation options
wf_method = 'spherical wavefront';
x_pixels = 750; % pixels per period for simulation
reptimes = 15; % repeat g1 so that large z aren't smeared
padsize = x_pixels*reptimes;

E_x = E_0;
E_spectrum = 1;
N = x_pixels*reptimes; 
lambda = lambda_from_E(E_x)';
k = 2*pi./lambda';

vis_map = zeros(length(z),1);
amp_map = zeros(length(z),1);
vis_map1 = zeros(length(z),1);
amp_map1 = zeros(length(z),1);
curves_mat = zeros(steps,length(z),1);
carpets_mat = zeros(N,length(z),1);
vis = zeros(1,length(z)); amp = zeros(1,length(z));
curves = zeros(steps,length(z));
g1_pattern = zeros(N,1);
wf1  = zeros(N,1);

tic
pixsize = p1/x_pixels;
x = (-(N/2):(N/2-1))*pixsize; % real space coordinates    
%% Create initial wavefront and grating    
wf1 = zeros(N,1);
[delta1,beta1] = get_refindex(m1, E_x);
[delta2,beta2] = get_refindex(m2, E_x);
tmp1 = E_spectrum*exp(-1i*(delta1-1i*beta1)*k*t1);
tmp2 = E_spectrum*exp(-1i*(delta2-1i*beta2)*k*t2);
temp = [tmp1*ones(1,round(x_pixels/2)) tmp2*ones(1,round(x_pixels/2))]';
temp = repmat(temp,reptimes,1);
g1_pattern = temp;
% rotate grating
t_pix = t1/pixsize;
unrot_g1 = zeros(t_pix,x_pixels*reptimes);
for i = 1:x_pixels:size(unrot_g1,2)
    unrot_g1(:,i:i+round(x_pixels/2)) = t1/(t_pix);
end
unrot_g1 = repmat(unrot_g1,1,2);
unrot_g1 = padarray(unrot_g1,[round(size(unrot_g1,1)/10) 0]);
rot_g1 = imrotate(unrot_g1,theta,'nearest','crop');
figure, imagesc(unrot_g1)
figure, imagesc(rot_g1)
figure, plot(sum(unrot_g1,1)), hold on, plot(sum(rot_g1,1))
clear unrot_g1
rot_g1 = sum(rot_g1,1);
rot_g1 = rot_g1(round(end/4):round(3*end/4)-1)';
temp = E_spectrum*exp(-1i*(delta1-1i*beta1)*k*rot_g1);


wf1  = temp.*exp( 1i.*k.*sqrt(d_sg1.^2+x.^2) )'./sqrt(d_sg1.^2+x.^2)';
wf1 = padarray(wf1,[padsize 0]);

%% Propagate wavefront
I_full = wf_propagate(wf1,z,E_x,pixsize,source_size,d_sg1, 'method','spherical wavefront');
I_full = I_full(padsize:(end-padsize-1),:);
%% phase stepping

for dis = 1:length(z)
    M = (d_sg1+z(dis))/d_sg1;
    p2 = p1*M/2;
    [stepping_curve]=phase_stepping(I_full(:,dis),pixsize,detector_pixsize,p2,0.5,steps,periods);
    ft_sc = fft(stepping_curve);
    vis(dis) = 2*abs(ft_sc(1+periods))./abs(ft_sc(1));
    amp(dis) = mean(stepping_curve);
    curves(:,dis) = stepping_curve;
end
toc

figure, imagesc(I_full)
figure, plot(z,vis)
figure, imagesc(curves)


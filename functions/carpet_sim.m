function [vis,amp,carpet] = carpet_sim(z,d_sg1,p1,dc,E_x,E_spectrum,m1,t1)
%CARPET_SIM Simulates talbot carpet from specified parameters
%   Inputs:
%           z: talbot distances [m]
%           d_sg1: source to g1 distance [m]
%           p1: period of first grating [m]
%           dc: duty cycle [0-1] (note: if DC = 0.4, this means m2 is wider than m1)
%           E_x: energies for energy spectrum [eV]
%           E_spectrum: intensities of source at energies given by E_x [1]
%           m1: first material of first grating (i.e. 'Au')
%           t1: thickness of first material of first grating [m]
%   Outputs:
%           vis: visibility for each propagation distance z
%           carpet: intensity pattern at each distance z

%% Hardcoded stuff
x_pixels = 50;
% calculate required reptimes
M = (d_sg1+z)/d_sg1;
reptimes = ceil(max(M(:)))+1;
%reptimes = 15;
t2 = 0;
m2 = 'Au';
padsize = round(x_pixels*reptimes/2);
source_size = 1.5e-06;
%steps = 15; % phase steps
periods = 2;
% % hardcoded 'spherical wavefront' option for initial wavefront
% % hardcoded 'spherical wavefront' option for propagation method
% % hardcoded 'direct' option for visibility calculation
% % using wf_propagate_hardcoded.m, which has certain options hardcoded
%% Frequently used values
N = x_pixels*reptimes; 
lambda = lambda_from_E(E_x)';
k = 2*pi./lambda';
pixsize = p1/x_pixels;
x = (-(N/2):(N/2-1))*pixsize; % real space coordinates    
%% initialize
vis = zeros(1,length(z)); 
amp = zeros(1,length(z));
%% Create initial wavefront and grating    
wf0 = exp( 1i.*k'.*sqrt(d_sg1.^2+x.^2) )'./sqrt(d_sg1.^2+x.^2)'; % before g1
g1_pattern = build_g1(m1,m2,t1,t2,dc,E_x,E_spectrum,x_pixels,reptimes);
wf1 = wf0.*g1_pattern; % after g1, using projection approximation
wf1 = padarray(wf1,[padsize 0]); % pad for simulation to avoid edge reflection
%% Propagate wavefront
I_full = wf_propagate_hardcoded(wf1,z,E_x,pixsize,source_size,d_sg1);
I_full = I_full(padsize:(end-padsize-1),:);
%% Visibility calculation
    for dis = 1:length(z)
        psize = M(dis)*x_pixels;
        curve = I_full(round(N/2-psize/2):round(N/2+psize/2),dis);
        ft_curve = fft(curve);
        vis(dis) = 2*abs(ft_curve(1+periods))./abs(ft_curve(1));
        amp(dis) = mean(curve);
        %curves(:,dis) = interp1(1:length(curve),curve,1:length(curve)/steps:length(curve))';
    end
carpet = I_full;

end


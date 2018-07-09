function [vis,carpet,time] = carpet_sim(z,d_sg1,p1,E_x,E_spectrum,m1,t1)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here


tic
%% Hardcoded stuff
x_pixels = 50;
reptimes = 5;
t2 = 0;
m2 = 'Au';
padsize = round(x_pixels*reptimes/2);
source_size = 2e-06;
steps = 15; % phase steps
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
%amp = zeros(1,length(z));
%curves = zeros(steps,length(z));
g1_pattern = zeros(N,length(E_x));
%% Create initial wavefront and grating    
wf1 = zeros(N,length(E_x));
for e = 1:length(E_x)
    [delta1,beta1] = get_refindex(m1, E_x(e));
    [delta2,beta2] = get_refindex(m2, E_x(e));
    tmp1 = E_spectrum(e)*exp(-1i*(delta1-1i*beta1)*k(e)*t1);
    tmp2 = E_spectrum(e)*exp(-1i*(delta2-1i*beta2)*k(e)*t2);
    temp = [tmp1*ones(1,round(x_pixels/2)) tmp2*ones(1,round(x_pixels/2))]';
    temp = repmat(temp,reptimes,1);
    g1_pattern(:,e) = temp;
    % % hardcoded 'spherical wavefront' option
        wf1(:,e)  = temp.*exp( 1i.*k(e).*sqrt(d_sg1.^2+x.^2) )'./sqrt(d_sg1.^2+x.^2)';
end
wf1 = padarray(wf1,[padsize 0]);
%% Propagate wavefront
I_full = wf_propagate_hardcoded(wf1,z,E_x,pixsize,source_size,d_sg1);
I_full = I_full(padsize:(end-padsize-1),:);
%% Visibility calculation
M = (d_sg1+z)/d_sg1;
    for dis = 1:length(z)
        psize = M(dis)*x_pixels;
        curve = I_full(round(N/2-psize/2):round(N/2+psize/2),dis);
        ft_curve = fft(curve);
        vis(dis) = 2*abs(ft_curve(1+periods))./abs(ft_curve(1));
        %amp(dis) = mean(curve);
        %curves(:,dis) = interp1(1:length(curve),curve,1:length(curve)/steps:length(curve))';
    end
carpet = I_full;
time = toc

end


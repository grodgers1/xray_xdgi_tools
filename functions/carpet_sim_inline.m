function [vis,amp,carpet] = carpet_sim_inline(z,d_sg1,p1,dc,E_x,E_spectrum,m1,t1)
%CARPET_SIM Simulates talbot carpet from specified parameters. This version
%has everything inline, to speed it up
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
lambda = (1.239842*10^(-6)./E_x)';
k = 2*pi./lambda';
pixsize = p1/x_pixels;
x = (-(N/2):(N/2-1))*pixsize; % real space coordinates    
%% initialize
vis = zeros(1,length(z)); 
amp = zeros(1,length(z));
%% Create initial wavefront and grating    
wf0 = exp( 1i.*k'.*sqrt(d_sg1.^2+x.^2) )'./sqrt(d_sg1.^2+x.^2)'; % before g1
% build g1 (inline version of build_g1.m)
sb1 = ceil(x_pixels*dc); % size of grating bar material 1 (pixels)
sb2 = x_pixels-sb1; % size of grating bar material 2 (pixels)
g1_pattern = zeros(N,length(E_x)); % complex transmission function
for e = 1:length(E_x)
    [delta1,beta1] = get_refindex(m1, E_x(e));
    [delta2,beta2] = get_refindex(m2, E_x(e));
    tmp1 = E_spectrum(e)*exp(-1i*(delta1-1i*beta1)*k(e)*t1);
    tmp2 = E_spectrum(e)*exp(-1i*(delta2-1i*beta2)*k(e)*t2);
    temp = [tmp1*ones(1,sb1) tmp2*ones(1,sb2)]';
    temp = repmat(temp,reptimes,1);
    g1_pattern(:,e) = temp;
end
% use projection approximation to get wf after g1
wf1 = wf0.*g1_pattern; % after g1
wf1 = padarray(wf1,[padsize 0]); % pad for simulation to avoid edge reflection
%% Propagate wavefront
% inline version of wf_propagate_hardcoded.m:
N2 = size(wf1,1); % length of wavefront
x2 = (-(N2/2):(N2/2-1))*pixsize; % real space coordinates
u  = (-(N2/2):(N2/2-1))./(N2*pixsize); % fourier space coordinates
f_wf1  = fftshift(fft(ifftshift(wf1,1),[],1),1); % ft wavefront (for convolution)
I = zeros(length(x2),length(E_x),length(z)); % initialize
for dis = 1:length(z)
    f_prop = exp(-1i*pi*lambda.*z(dis).*u.^2); % ft Fresnel propagator
    f_prop = f_prop';
    norm_factor = exp(1i*k.*z(dis))./(2*d_sg1); % normalization
    norm_factor = norm_factor./abs(norm_factor);
    p_wf = norm_factor.*ifftshift(ifft(f_wf1.*f_prop,[],1),1); % convolute
    p_I = abs(p_wf).^2; % propagated intensity
        % source size considerations
        w = z(dis)*(source_size/d_sg1); % Consider magnification of source
        prof_source = exp(-(1/2)*(x2.^2)/ w^2)'; % gaussian source
        prof_source = prof_source./sum(prof_source); % normalize
        ft_source = fftshift(fft(prof_source,[],1),1); % ft source for convolution
        ft_I = fftshift(fft(p_I,[],1),1); % ft wavefront intensity for convolution
        temp = abs(ifftshift(ifft(ft_I.*ft_source,[],1),1)); % convolution
        I(:,:,dis) = temp;

end
I_full = squeeze(sum(I,2)); % incoherently sum different energy contributions
%
I_full = I_full(padsize:(end-padsize-1),:); % crop out padding
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


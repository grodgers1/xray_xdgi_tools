function [I_z] = wf_propagate_hardcoded(wf0,z,E,dx,width,d_source)
% wf_propagate Propagates an initial wavefront to positions given by z
%    Uses the Fresnel equation and paraxial approximation
% Inputs:
%        wf0: initial wavefront 
%            (complex array with size N-by-length(E))
%        E: Energies
%        dx: [m] pixel size (sampling rate) of wf0 
%        width: [m] source size (default is fwhm) 
%        d_source: [m] distance of source to wf0.
% Outputs:
%       I_z: Intensity of propagated wavefront

lambda = lambda_from_E(E)';
k = 2*pi./lambda';

N = size(wf0,1); % length of wavefront
x = (-(N/2):(N/2-1))*dx; % real space coordinates
u  = (-(N/2):(N/2-1))./(N*dx); % fourier space coordinates

f_wf0  = fftshift(fft(ifftshift(wf0,1),[],1),1); % ft wavefront (for convolution)
I = zeros(length(x),length(E),length(z)); % initialize
for dis = 1:length(z)
    f_prop = exp(-1i*pi*lambda.*z(dis).*u.^2); % ft Fresnel propagator
    f_prop = f_prop';
    norm_factor = exp(1i*k.*z(dis))./(2*d_source); % normalization
    norm_factor = norm_factor./abs(norm_factor);
    p_wf = norm_factor.*ifftshift(ifft(f_wf0.*f_prop,[],1),1); % convolute
    p_I = abs(p_wf).^2; % propagated intensity

        w = z(dis)*(width/d_source); % Consider magnification of source
        prof_source = exp(-(1/2)*(x.^2)/ w^2)'; % gaussian source
        prof_source = prof_source./sum(prof_source); % normalize
        ft_source = fftshift(fft(prof_source,[],1),1); % ft source for convolution
        ft_I = fftshift(fft(p_I,[],1),1); % ft wavefront intensity for convolution
        temp = abs(ifftshift(ifft(ft_I.*ft_source,[],1),1)); % convolution
        I(:,:,dis) = temp;

end
I_z = squeeze(sum(I,2)); % incoherently sum different energy contributions

end


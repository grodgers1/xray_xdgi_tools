function [ I_z, WF_z ] = fresnel_propagator( WF0, pixsize, z, lambda, method )
%fresnel_propagator Gives the wavefront a distance z downstream
%   I0: initial intensity profile of wavefront
%   p0: initial phase profile of wavefront
%   resolution: pixel size of initial (and final) profiles [m]
%   z: distance downstream to propagate [m]
%   lamda: wavelength [m]
%   method: how to perform convolution with fresnel propagator?
%   1: real space prop and conv with fft. 2: k space prop and conv with fft
%   3: real space prop and conv with matlab's conv() function (not yet implemented)
k=2*pi/lambda; % get wavenumber
sx = length(WF0); % size of input wavefunction
x = ((0:sx-1)-(sx-1)/2)*pixsize; % position vector
norm_factor = mean(mean(real(WF0))); % for normalizing intensity of carpet
I_z = zeros(length(z),length(WF0)); % initialize Talbot carpet (intensity profile)
WF_z = zeros(length(z),length(WF0)); % initialize wavefront at position z

% Method 1: real space propagator, convoluted by fft then multiplication
% with fft of wavefunction
fft_WF = fftshift(fft(WF0)); % fourier transfer for convolution
if method == 1
    for ii = 1:length(z)
        fprintf(['Getting wavefront for step ' num2str(ii) ' of ' num2str(length(z)) '\n'])
        prop = exp(1i*k*z(ii)).*exp(-1i*k*(x.^2)/(2*z(ii)))/(2*pi*z(ii));  % fresnel propagator
        fft_prop = fftshift(fft(prop)); % fourier transfer propagator for convolution
        WF_propagated = ifft(ifftshift(fft_WF.*fft_prop));  % convolute propagator with initial wavefunction
        I_z(ii,:) = (abs(WF_propagated)').^2; % get intensity of wavefront at position z
        div_factor = mean(I_z(ii,:));
        I_z(ii,:) = norm_factor * I_z(ii,:)/div_factor;
        WF_z(ii,:) = WF_propagated; % Wavefunction at position z
    end
% Method 2: fourier space propagator, convoluted by multiplication with fft
% of wavefunction
elseif method == 2
    k_x = 1./x; % reciporical space conjugate to x
    for ii = 1:length(z)
        fprintf(['Getting wavefront for step ' num2str(ii) ' of ' num2str(length(z)) '\n'])
        fprop = fftshift(exp(-1i*pi*lambda*z(ii)*(k_x.^2)));  % fresnel propagator in fourier space
        WF_propagated = ifft(ifftshift(fft_WF.*fprop));  % convolute propagator with initial wavefunction
        I_z(ii,:) = (abs(WF_propagated)').^2; % get intensity of wavefront at position z
        div_factor = mean(I_z(ii,:));
        I_z(ii,:) = norm_factor * I_z(ii,:)/div_factor;
        WF_z(ii,:) = WF_propagated; % Wavefunction at position z
    end
    fprintf('This method may have problems right now')
elseif method == 3
    % % Method 3: convoltion using matlab's conv() function
    % WF_0 = I0.*exp(1i*p0); % initial complex wavefunction
    % for ii = 1:length(z)
    %     fprintf(['Getting wavefront for step ' num2str(ii) ' of ' num2str(length(z)) '\n'])
    %     prop = exp(1i*k*z(ii)).*exp(-1i*k*(x.^2)/(2*z(ii)))/(2*pi*z(ii));  % fresnel propagator in real space
    %     WF_propagated = conv(WF0,prop); % convolute propagator with initial wavefunction
    %     WF_propagated = WF_propagated(:); % cut to proper size 
    %     I_z(ii,:) = (abs(WF_propagated)').^2; % get intensity of wavefront at position z
    %     div_factor = mean(I_z(ii,:));
    %     I_z(ii,:) = norm_factor * I_z(ii,:)/div_factor;
    %     WF_z(ii,:) = WF_propagated; % Wavefunction at position z
    % end
    fprintf('This method is not yet implemented')
end
end


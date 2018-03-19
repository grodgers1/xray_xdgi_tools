%% Builds talbot carpet with options for the following:
% - Any intensity and phase wavefront from g1 (done)
% - Automatic g1 wavefront from materials descritpion (need to do)
% - Gaussian source (need to do)
% - Cone beam (need to do)
% - Polychromatic spectrum (need to do)
% - Phase stepping over pattern (done)

%% Geometry and setup parameters:

% g1 details
p1 = 7e-06; % [m] period of the above mentioned pattern
m1 = 'Si'; % First material (e.g. 'Si')
m2 = 'Au'; % Second material (e.g. 'Au')
t1 = 6e-06; % thickness of first material [m]
t2 = 6e-06; % thickness of second material [m]
z = (0.01:0.01:4)*(7e-06)^2/(4*lambda);

% p2 and detector details
p2 = p1/2; % period of 2nd grating [m] (p1/2 for pi shift)
z_detector = 150; % position of detector [index of z array]
steps = 15; % How many phase steps?
periods = 1; % How many periods to phase step over?

% Source details
E_spectrum = (20:30)*1e+03; % Energies for above intensity spectrum [eV]
I_spectrum = ones(size(E_spectrum)); % Intensity at energy described by E_spectrum
% source_size = 2e-06; % FWHM of source [m]
% d_sg1 = 0.3; % source to g1 distance [m]

%% Computing and visualizing options:
x_pixels = 2000; % pixels per period for simulation
reptimes = 20; % repeat g1 so that large z aren't smeared
method = 1; % see fresnel_propagator.m for info
periods_to_plot = 3;

% %% Find magnification
% M = (d_sg1+z(z_detector))/d_sg1; % magnification factor

%% Building initial wavefunction (no more inputs needed)
tmp = round(x_pixels/2);
I_pattern = zeros(2*tmp, length(E_spectrum));
for e = 1:length(E_spectrum)
    lambda = lambda_from_E(E_spectrum(e));
    k = 2*pi./lambda;
    delta0 = find_delta(m1, lambda, 'wavelength');
    delta1 = find_delta(m2, lambda, 'wavelength');
    beta0 = find_beta(m1, lambda, 'wavelength');
    beta1 = find_beta(m2, lambda, 'wavelength');
    I_pattern(:,e) = [exp(-2*k*beta0*t1)*ones(tmp,1)' exp(-2*k*beta1*t2)*ones(tmp,1)'];
    p_pattern(:,e) = [-k*t1*delta0*ones(tmp,1)' -k*t2*delta1*ones(tmp,1)'];
end

WF0 = I_pattern.*exp(1i*p_pattern);
% enlarge wavefront for large grating
WF0 = repmat(WF0,reptimes);
WF0 = WF0(1:length(E_spectrum),:);

% %% Convolute wavefunction with source
% sz_WF0 = length(WF0);
% x = -((sz_WF0-1)/2):((sz_WF0-1)/2);
% source = exp(-(x./source_size).^2); % create gaussian source

%% Running propagator to get wavefunction at a distance z downstream:
pixsize = p1/x_pixels;
for e = 1:length(E_spectrum)
    [I_e_z(:,e),WF_e_z(:,e)] = fresnel_propagator(WF0(:,e), pixsize, z, lambda_from_E(E_spectrum(e)), method);
end

WF_z = sum(WF_e_z,2); % Sum up wavefunction from various wavelength contributions
I_z = abs(WF_z); % Find intensity of full wavefunction

%% Plotting results
temp = periods_to_plot*x_pixels;
temp2 = round(size(I_z,2)/2)-round(temp/2);
figure, imagesc(I_z(:,temp2:(temp2+temp))'), colormap gray

%% Phase stepping:
[stepping_curve,phase,ac,vis]=phase_stepping(I,pixsize,detector_pixsize,p2,steps,periods);
figure, plot(stepping_curve), xlabel('phase step'), ylabel('intensity')











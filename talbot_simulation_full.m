%% Builds talbot carpet with options for the following:
% - Any intensity and phase wavefront from g1 (done)
% - Automatic g1 wavefront from materials descritpion (need to do)
% - Gaussian source (need to do)
% - Cone beam (need to do)
% - Polychromatic spectrum (need to do)
% - Phase stepping over pattern (done)

%% Geometry and setup parameters:
lambda = lambda_from_E(20); % [m] (use lambda_from_E(20) to input E in keV)
I_pattern = [1 1];
p_pattern = [pi 0];
p1 = 7e-06; % [m] period of the above mentioned pattern
% dc_g1 = 0.5; % duty cycle of first grating 
x_pixels = 2000; % pixels per period for simulation
repg1 = 20; % repeat g1 so that large z aren't smeared
z = (0.01:0.01:4)*(7e-06)^2/(4*lambda);

p2 = p1/2; % [m] period of 2nd grating (p1/2 for pi shift)
z_detector = 150; % position of detector [index of z array]
steps = 15;
periods = 1;
detector_pixsize = 10e-06;

% t0 = 6e-06;
% t1 = 6e-06;
% m0 = 'Si';
% m1 = 'Au';
% delta0 = find_delta(m0, lambda, 'wavelength');
% delta1 = find_delta(m1, lambda, 'wavelength');
% beta0 = find_beta(m0, lambda, 'wavelength');
% beta1 = find_beta(m1, lambda, 'wavelength');
% 
% I_pattern = [exp(-mu0*t0) exp(-mu1*t1)];
% p_pattern = [t0*delta0 t1*delta1];

%% Computing and visualizing options:
method = 1; % see fresnel_propagator.m for info

periods_to_plot = 3;

%% Setting up wave fronts:
tmp = round(x_pixels/2);
I0 = [I_pattern(1)*ones(tmp,1)' I_pattern(2)*ones(tmp,1)'];
I0 = repmat(I0,repg1);
I0 = I0(1,:);
p0 = [p_pattern(1)*ones(tmp,1)' p_pattern(2)*ones(tmp,1)'];
p0 = repmat(p0,repg1);
p0 = p0(1,:);

%% Running propagator:
pixsize = p1/x_pixels;
[I_z,~] = fresnel_propagator(I0, p0, pixsize, z, lambda, method);
temp = periods_to_plot*x_pixels;
temp2 = round(size(I_z,2)/2)-round(temp/2);
figure, imagesc(I_z(:,temp2:(temp2+temp))'), colormap gray

%% Phase stepping:
[stepping_curve,phase,ac,vis]=phase_stepping(I,pixsize,detector_pixsize,p2,steps,periods);
figure, plot(stepping_curve), xlabel('phase step'), ylabel('intensity')


% make talbot carpet:
pattern = [zeros(1000,1)' pi*ones(1000,1)'];
period = length(pattern);
p0 = repmat(pattern,20);
p0 = p0(1,:);
I0=ones(size(p0));
p1 = 7e-06; % grating 1 period [m]
pixsize = (p1/period);
lambda = lambda_from_E(20);
z = (0.01:0.01:4)*(7e-06)^2/(4*lambda);
tic
[I_z_phase,~] = fresnel_propagator(I0, p0, pixsize, z, lambda, 1);
toc
carpet_phase = [I0' I_z_phase'];
temp = round(size(carpet_phase,1)/2);
figure, imagesc(carpet_phase(temp:(temp+3*period),:)), colormap gray

% step over carpet
I = I_z_phase(50,:);
steps = 15;
periods = 1;
p2 = (7e-06)/2; % p2 = p1/2 for pi shift, p2 = p1 for abs or pi/2 shift
detector_pixsize = 10e-06;
[stepping_curve,phase,ac,vis]=phase_stepping(I,pixsize,detector_pixsize,p2,steps,periods);
figure, plot(stepping_curve), xlabel('phase step'), ylabel('intensity')

% Alternative by convolution with grating 2 (Harti et al. 2017 Optics Express)
% construct the absorption grating
halfperiod_g2 = round(period/4);
g2_pattern = [zeros(halfperiod_g2,1)' 1*ones(halfperiod_g2,1)'];
g2_period = length(g2_pattern);
g2 = g2_pattern;
% g2 = repmat(g2_pattern,20);
% g2 = g2(1,:);
I_detector = conv(I, g2, 'same');



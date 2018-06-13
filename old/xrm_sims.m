%% Simulations for xrm abstract

E = 30000;
lambda = lambda_from_E(E);

p1 = 7e-06; % [m] period of the above mentioned pattern
m1 = 'Au'; % First material (e.g. 'Si')
% m2 = 'Au'; % Second material (e.g. 'Au')
t1 = 6e-06; % thickness of first material [m]
t2 = 6e-06; % thickness of second material [m]
z = (0.01:0.01:3)*(7e-06)^2/(4*lambda);


x_pixels = 2000; % pixels per period for simulation
reptimes = 20; % repeat g1 so that large z aren't smeared
method = 1; % see fresnel_propagator.m for info
periods_to_plot = 3;



tmp = round(x_pixels/2);
lambda = lambda_from_E(E);
k = 2*pi./lambda;
[delta0, beta0] = get_refindex(m1, lambda);

I_pattern = [exp(-2*k*beta0*t1)*ones(tmp,1)' ones(tmp,1)'];
p_pattern = [k*t1*delta0*ones(tmp,1)' zeros(tmp,1)'];


WF0 = I_pattern.*exp(1i*p_pattern);
% figure, plot(real(WF0))
% hold on, plot(imag(WF0))

% enlarge wavefront for large grating
WF0 = repmat(WF0,reptimes);
WF0 = WF0(1,:);

pixsize = p1/x_pixels;

[I_z,WF_z] = fresnel_propagator(WF0, pixsize, z, lambda, method);

temp = periods_to_plot*x_pixels;
temp2 = round(size(I_z,2)/2)-round(temp/2);
figure, imagesc(I_z(:,temp2:(temp2+temp))'), colormap gray

I_pattern2 = [ones(tmp,1)' ones(tmp,1)'];
p_pattern2 = [k*t1*delta0*ones(tmp,1)' zeros(tmp,1)'];

WF2 = I_pattern2.*exp(1i*p_pattern2);
% figure, plot(real(WF2))
% hold on, plot(imag(WF2))

% enlarge wavefront for large grating
WF2 = repmat(WF2,reptimes);
WF2 = WF2(1,:);

[I_z2,WF_z2] = fresnel_propagator(WF2, pixsize, z, lambda, method);

temp = periods_to_plot*x_pixels;
temp2 = round(size(I_z2,2)/2)-round(temp/2);
figure, imagesc(I_z2(:,temp2:(temp2+temp))'), colormap gray


I_pattern3 = [zeros(tmp,1)' ones(tmp,1)'];
p_pattern3 = [zeros(tmp,1)' zeros(tmp,1)'];

WF3 = I_pattern3.*exp(1i*p_pattern3);
% figure, plot(real(WF3))
% hold on, plot(imag(WF3))

% enlarge wavefront for large grating
WF3 = repmat(WF3,reptimes);
WF3 = WF3(1,:);

[I_z3,WF_z3] = fresnel_propagator(WF3, pixsize, z, lambda, method);

temp = periods_to_plot*x_pixels;
temp2 = round(size(I_z3,2)/2)-round(temp/2);
figure, imagesc(I_z3(:,temp2:(temp2+temp))'), colormap gray


figure, plot(I_z(240,temp2:(temp2+temp)))
hold on, plot(I_z2(240,temp2:(temp2+temp)))
hold on, plot(I_z3(240,temp2:(temp2+temp)))

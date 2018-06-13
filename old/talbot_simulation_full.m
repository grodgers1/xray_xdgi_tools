%% Builds talbot carpet with options for the following:
% - Any intensity and phase wavefront from g1 (done)
% - Automatic g1 wavefront from materials descritpion (done)
% - Gaussian source (need to do)
% - Cone beam (need to do)
% - Polychromatic spectrum (done)
% - Phase stepping over pattern (done)

addpath('/Users/Griffin/Documents/MATLAB/xray_xdgi_tools/tables')
addpath('/Users/Griffin/Documents/MATLAB/xray_xdgi_tools')

%% Geometry and setup parameters:

% g1 details
p1 = 7e-06; % [m] period of the above mentioned pattern
m1 = 'Au'; % First material (e.g. 'Si')
m2 = 'Au'; % Second material (e.g. 'Au')
t1 = 6e-06; % thickness of first material [m]
t2 = 0; % thickness of second material [m]

z = (0.01:0.01:4)*(7e-06)^2/(4*lambda_from_E(30000));

% p2 and detector details
p2 = p1/2; % period of 2nd grating [m] (p1/2 for pi shift)
z_detector = 150; % position of detector [index of z array]
steps = 15; % How many phase steps?
periods = 1; % How many periods to phase step over?

% Source details
% % Gaussian source:
% E_spectrum = (25:0.5:35)*1e+03; % Energies for above intensity spectrum [eV]
% I_spectrum = exp(-((E_spectrum-30000)/4000).^2); % Intensity at energy described by E_spectrum
% figure, plot(E_spectrum, I_spectrum)
% % Monochromatic source
E_spectrum = 30000;
I_spectrum = 1;

% Source size/location
source_size = 2e-06; % FWHM of source [m]
d_sg1 = 0.294; % source to g1 distance [m]

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
p_pattern = zeros(2*tmp, length(E_spectrum));
k = 2*pi./lambda_from_E(E_spectrum);

for e = 1:length(E_spectrum)
    [delta1,beta1] = get_refindex(m1, E_spectrum(e));
    [delta2,beta2] = get_refindex(m2, E_spectrum(e));
    I_pattern(:,e) = [I_spectrum(e)*exp(-2*k*beta1*t1)*ones(tmp,1)' I_spectrum(e)*exp(-2*k*beta2*t2)*ones(tmp,1)'];
    p_pattern(:,e) = [-k*t1*delta1*ones(tmp,1)' -k*t2*delta2*ones(tmp,1)'];
end

WF0 = I_pattern.*exp(1i*p_pattern);
% enlarge wavefront for large grating
WF0 = repmat(WF0,reptimes);
WF0 = WF0(:,1:length(E_spectrum));

if plots_on ==1
    figure, plot(real(WF0(1:x_pixels,:)),'r')
    hold on
    plot(imag(WF0(1:x_pixels,:)),'b')
    plot(abs(WF0(1:x_pixels,:)).^2,'k')
    legend({'Re{WF}','Im{WF}','Intensity'})
end

%% Running propagator to get wavefunction at a distance z downstream:
pixsize = p1/x_pixels;
I_e_z = zeros(length(WF0),length(z),length(E_spectrum));
% WF_e_z = zeros(length(WF0),length(z),length(E_spectrum));
for e = 1:length(E_spectrum)
    fprintf(['Energy step ' num2str(e) ' of ' num2str(length(E_spectrum)) '\n'])
    [I_e_z(:,:,e),~] = fresnel_propagator(WF0(:,e), pixsize, z, lambda_from_E(E_spectrum(e)), method);
end

I_z = sum(I_e_z,3);
%clear I_e_z

%% Convoluting with source spot size
% x = -((length(WF0)-1)/2):((length(WF0)-1)/2);
% source = zeros(length(WF0),length(z));
% for dist = 1:length(z)
% source(:,dist) = exp(-(x*d_sg1/(z(dist)*source_size)).^2);
% end
% 
% I_zf = ifft(fft(source).*fft(I_z));


%% Plotting results
temp = periods_to_plot*x_pixels;
temp2 = round(size(I_z,1)/2)-round(temp/2);
figure, imagesc(I_z(temp2:(temp2+temp),:)), colormap gray

        


%% Phase stepping:
detector_pixsize = 1;
steps = 10;
[stepping_curve,phase,ac,vis]=phase_stepping(I_z(:,400),pixsize,detector_pixsize,p2,steps,periods,1);
figure, plot(stepping_curve), xlabel('phase step'), ylabel('intensity')











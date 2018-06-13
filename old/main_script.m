%% Builds talbot carpet with options for the following:
% - Any intensity and phase wavefront from g1 (done)
% - Automatic g1 wavefront from materials descritpion (done)
% - Gaussian source (need to do)
% - Cone beam (need to do)
% - Polychromatic spectrum (done)
% - Phase stepping over pattern (done)
clear all
addpath('/Users/Griffin/Documents/MATLAB/xray_xdgi_tools/tables')
addpath('/Users/Griffin/Documents/MATLAB/xray_xdgi_tools')

%% Geometry and setup parameters:

% g1 details
p1 = 7e-06; % [m] period of the above mentioned pattern
m1 = 'Au'; % First material (e.g. 'Au')
m2 = 'Au'; % Second material (e.g. 'Au')
t1 = 6e-06; % thickness of first material [m]
t2 = 0; % thickness of second material [m]

%z = (0.01:0.01:4)*(7e-06)^2/(4*lambda_from_E(30000));

z = 0.005:0.005:1.2;
% % % recall d_t = 

% p2 and detector details
p2 = p1/2; % period of 2nd grating [m] (p1/2 for pi shift)

% Source details
source_spectrum = 1;

if source_spectrum == 0
% % Gaussian source:
E_middle = 30000;
E_width = 4000;
E_spectrum = (20:1:40)*1e+03; % Energies for above intensity spectrum [eV]
I_spectrum = exp(-((E_spectrum-E_middle)/E_width).^2); % Intensity at energy described by E_spectrum
figure, plot(E_spectrum, I_spectrum)
elseif source_spectrum == 1
% % Monochromatic source
E_spectrum = 30000;
I_spectrum = 1;
elseif source_spectrum == 2
% % Gaussian source:
E_middle = 35000;
E_width = 5000;
E_spectrum = (25:1:45)*1e+03; % Energies for above intensity spectrum [eV]
I_spectrum = exp(-((E_spectrum-E_middle)/E_width).^2); % Intensity at energy described by E_spectrum
figure, plot(E_spectrum, I_spectrum)
end

% Source size/location
source_size = 2e-06; % FWHM of source [m]
d_sg1 = 0.294; % source to g1 distance [m]

%% Computing and visualizing options:
plots_on = 1;
x_pixels = 2000; % pixels per period for simulation
reptimes = 20; % repeat g1 so that large z aren't smeared
method = 1; % see fresnel_propagator.m for info

pixsize = p1/x_pixels;

%% Building initial wavefunction
x = (-((x_pixels*reptimes-1)/2):((x_pixels*reptimes-1)/2))*pixsize;
k = 2*pi./lambda_from_E(E_spectrum);

tmp = round(x_pixels/2);
WF_pattern = zeros(2*tmp,length(E_spectrum));
for e = 1:length(E_spectrum)
    [delta1,beta1] = get_refindex(m1, E_spectrum(e));
    [delta2,beta2] = get_refindex(m2, E_spectrum(e));
    tmp1 = I_spectrum(e)*exp(-1i*(delta1-1i*beta1)*k(e)*t1);
    tmp2 = I_spectrum(e)*exp(-1i*(delta2-1i*beta2)*k(e)*t2);
    WF_pattern(:,e) = [tmp1*ones(1,tmp) tmp2*ones(1,tmp)]';
end
% Make wavefront bigger than 1 period
WF_pattern = repmat(WF_pattern,reptimes,1);
% % find wavefront after 1st grating using projection approximation
%WF_1 = WF_pattern.*WF0;
WF_1 = WF_pattern;

if plots_on ==1
    figure, plot(real(WF_1),'r')
    hold on
    plot(imag(WF_1),'b')
    plot(abs(WF_1).^2,'k')
    legend({'Re{WF}','Im{WF}','Intensity'})
end

%% Running propagator to get wavefront at a distance z downstream:
I_e_z = zeros(length(WF_1),length(z),length(E_spectrum));
% WF_e_z = zeros(length(WF0),length(z),length(E_spectrum));
for e = 1:length(E_spectrum)
    fprintf(['Energy step ' num2str(e) ' of ' num2str(length(E_spectrum)) '\n'])
    [I_e_z(:,:,e),~] = fresnel_propagator(WF_1(:,e), pixsize, z, lambda_from_E(E_spectrum(e)), method);
end

I_z = sum(I_e_z,3);
clear I_e_z

%% Convoluting with source spot size
% % convolution with wavefront
if plots_on == 1
    figure, hold on,
    xlabel('position [m]'), ylabel('source profile [intensity]')
end
ft_I_z = fft(I_z);
for pos = 1:length(z)
    %tic
    % % find projected source prof (considering magnification)
    source_prof = exp(-x.^2/(2*(source_size*(z(pos)/d_sg1))^2));
    I_z(:,pos) = ifft(ft_I_z(:,pos).*fft(source_prof)');
    %%%%%% Need to add padding for fft
    if mod(pos,25) == 0 && plots_on == 1
        plot(x,source_prof)
    end
    %toc
end

%% Plotting results
periods_to_plot = 5;
temp = periods_to_plot*x_pixels;
temp2 = round(size(I_z,1)/2)-round(temp/2);
figure, imagesc(I_z(temp2:(temp2+temp),:)), colormap gray
hold on
xlabel('inter-grating distance [m]')
ylabel('x position [m]')

figure, plot(I_z(temp2:(temp2+temp),round(end/2)))


%% Phase stepping:
detector_pixsize = 1;
steps = 15; % How many phase steps?
periods = 2; % How many periods to phase step over?

[stepping_curve,~,~,~]=phase_stepping(I_z(:,round(end/2)),pixsize,detector_pixsize,p2,0.5,steps,periods);
if plots_on == 1
figure, plot(stepping_curve,'o-')
xlabel('phase step'), ylabel('intensity')
end

%%
vis = zeros(1,length(z));
figure, hold on, xlabel('phase step'), ylabel('intensity')
for pos = 1:length(z)
    [stepping_curve,~,~,~]=phase_stepping(I_z(:,pos),pixsize,detector_pixsize,p2,0.5,steps,periods);
    vis(pos) = (max(stepping_curve)-min(stepping_curve))/(max(stepping_curve)+min(stepping_curve));
    if mod(pos,10) == 0
        plot(stepping_curve, 'o-')
    end
end

if plots_on == 1
    figure, plot(z,vis)
    xlabel('g2 position [m]')
    ylabel('visibility')
    ylim([0 1])
end




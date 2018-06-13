clear all

addpath('/Users/Griffin/Documents/MATLAB/xray_xdgi_tools/tables')
addpath('/Users/Griffin/Documents/MATLAB/xray_xdgi_tools')


%% Geometry and setup parameters:

% g1 details
m1 = 'Au'; % First material (e.g. 'Au')
m2 = 'Au'; % Second material (e.g. 'Au')
t1 = 6e-06; % thickness of first material [m]
t2 = 0; % thickness of second material [m]

z = 0.005:0.005:0.6;

% Source details
source_spectrum = 'monochromatic';

if source_spectrum == 'gaussian'
% % Gaussian source:
E_middle = 30000;
E_width = 4000;
E_spectrum = (20:1:40)*1e+03; % Energies for above intensity spectrum [eV]
I_spectrum = exp(-((E_spectrum-E_middle)/E_width).^2); % Intensity at energy described by E_spectrum
figure, plot(E_spectrum, I_spectrum)
elseif source_spectrum == 'monochromatic'
% % Monochromatic source
E_spectrum = 30000;
I_spectrum = 1;
elseif source_spectrum == 'tube'
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

%% Run loop over p1
p1_range = (2.4:0.3:8)*1e-06;
vis_mat = zeros(length(z),length(p1_range));
for p = 1:length(p1_range)
    fprintf(['Period ' num2str(p) ' of ' num2str(length(p1_range)) '\n'])
    tic
    p1 = p1_range(p); % [m] period of the above mentioned pattern
    p2 = p1/2;
    % % Building initial wavefunction
    tmp = round(x_pixels/2);
    k = 2*pi./lambda_from_E(E_spectrum);
    WF_pattern = zeros(2*tmp,length(E_spectrum));
    for e = 1:length(E_spectrum)
        [delta1,beta1] = get_refindex(m1, E_spectrum(e));
        [delta2,beta2] = get_refindex(m2, E_spectrum(e));
        tmp1 = I_spectrum(e)*exp(1i*(1-delta1+1i*beta1)*k(e)*t1);
        tmp2 = I_spectrum(e)*exp(1i*(1-delta2+1i*beta2)*k(e)*t2);
        WF_pattern(:,e) = [tmp1*ones(1,tmp) tmp2*ones(1,tmp)]';
    end
    % Make wavefront bigger than 1 period
    WF0 = repmat(WF_pattern,reptimes,1);

    % % Running propagator to get wavefront at a distance z downstream:
    pixsize = p1/x_pixels;
    I_e_z = zeros(length(WF0),length(z),length(E_spectrum));
    for e = 1:length(E_spectrum)
        fprintf(['Energy step ' num2str(e) ' of ' num2str(length(E_spectrum)) '\n'])
        [I_e_z(:,:,e),~] = fresnel_propagator(WF0(:,e), pixsize, z, lambda_from_E(E_spectrum(e)), method);
    end
    I_z = sum(I_e_z,3); clear I_e_z
    % % % axis for creating source profile
    x = (-((length(WF0)-1)/2):((length(WF0)-1)/2))*pixsize;
    % % convolution with wavefront
    ft_I_z = fft(I_z);
    for pos = 1:length(z)
        % % find projected source prof (considering magnification)
        source_prof = exp(-x.^2/(2*(source_size*(z(pos)/d_sg1))^2));
        I_z(:,pos) = ifft(ft_I_z(:,pos).*fft(source_prof)');
    end

    % % Show talbot carpet
%     periods_to_plot = 3;
%     temp = periods_to_plot*x_pixels;
%     temp2 = round(size(I_z,1)/2)-round(temp/2);
%     figure, imagesc(I_z(temp2:(temp2+temp),:)), colormap gray
%     hold on
%     xlabel('inter-grating distance [m]')
%     ylabel('x position [m]')
%     xticks(20:20:120)
%     xticklabels(round(z(20:20:120),2))
%     yticks(1000:1000:6000)
%     yticklabels((1000:1000:6000)*pixsize)
%     
    % % phase stepping options
    steps = 15; % How many phase steps?
    periods = 2; % How many periods to phase step over?
    % % Making vis vs. d plot
    vis = zeros(1,length(z));
    for pos = 1:length(z)
        [stepping_curve,~,~,~]=phase_stepping(I_z(:,pos),pixsize,1,p2,0.5,steps,periods);
        vis(pos) = (max(stepping_curve)-min(stepping_curve))/(max(stepping_curve)+min(stepping_curve));
    end
    vis_mat(:,p) = vis;
    toc
end

figure, plot(z, vis_mat(:,1:2:end))
xlabel('inter-grating distance [m]')
ylabel('visibility')


figure, imagesc(vis_mat')
xticks(20:20:120)
xticklabels(round(z(20:20:120),2))
yticks(1:5:length(p1_range))
yticklabels(p1_range(1:5:length(p1_range)))
xlabel('inter-grating distance [m]')
ylabel('period of g1 [m]')

% % find d value of each p1 choice
d_vec = zeros(1,length(p1_range));
for p = 1:length(p1_range)
    [pks, loc] = findpeaks();
    
    
    
    
end


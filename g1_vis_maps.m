
clear all
addpath('./tables/')

%% Settings
% % Define first grating
p1_range = (1:0.25:10)*1e-06; % [m] period of the above mentioned pattern
z = 0.0025:0.0025:0.3; % propagation distances [m]
steps = 15; % phase steps
m1 = 'Au'; % First material (e.g. 'Au')
m2 = 'Au'; % Second material (e.g. 'Au')
t1 = 6e-06; % thickness of first material [m]
t2 = 0; % thickness of second material [m]

% % Define source
% Size
source_size = 2e-06; % FWHM of source [m]
% Location
d_sg1 = 0.294; % source to g1 distance [m]
% Spectrum
source_spectrum = 0;
E_0 = 30000;
sig_E = 4000;
n = 25;
E_min = 15000;
E_max = 45000;
if source_spectrum == 0
% Gaussian source:
[E_spectrum,E_x] = EspectrumGauss(E_0, sig_E, n,E_min,E_max);
elseif source_spectrum == 1
% Monochromatic source
E_x = E_0;
E_spectrum = 1;
end
figure, plot(E_x, E_spectrum, '.-'), xlabel('Energy [eV]')
ylabel('Intensity'), title('Energy Spectrum')

% % Simulation options
spherical_wf = 0; % 1 = spherical wf, 0 = magnification. (if d_sg1 \approx z(end), you get artifacts from spherical (not sure why yet))
x_pixels = 750; % pixels per period for simulation
reptimes = 15; % repeat g1 so that large z aren't smeared

lambda = lambda_from_E(E_x)';
k = 2*pi./lambda';

vis_map = zeros(length(z),length(p1_range));
amp_map = zeros(length(z),length(p1_range));
vis_map1 = zeros(length(z),length(p1_range));
amp_map1 = zeros(length(z),length(p1_range));
curves_mat = zeros(steps,length(z),length(p1_range));
for p = 1:length(p1_range)
    
    tic
    p1 = p1_range(p);

    pixsize = p1/x_pixels;
    N = x_pixels*reptimes;
    % real space coordinates
    x = (-(N/2):(N/2-1))*pixsize;
    % fourier space coordinates
    u  = (-(N/2):(N/2-1))./(N*pixsize);

    %% Create initial wavefront and grating
    g1_pattern = zeros(length(x),length(E_x));
    wf1  = zeros(length(x),length(E_x));
    for e = 1:length(E_x)
        [delta1,beta1] = get_refindex(m1, E_x(e));
        [delta2,beta2] = get_refindex(m2, E_x(e));
        tmp1 = E_spectrum(e)*exp(-1i*(delta1-1i*beta1)*k(e)*t1);
        tmp2 = E_spectrum(e)*exp(-1i*(delta2-1i*beta2)*k(e)*t2);
        temp = [tmp1*ones(1,round(x_pixels/2)) tmp2*ones(1,round(x_pixels/2))]';
        temp = repmat(temp,reptimes,1);
        g1_pattern(:,e) = temp;
        if spherical_wf == 1
            wf1(:,e)  = temp.*(exp( 1i.*k(e).*sqrt(d_sg1.^2+x.^2) ))';
        elseif spherical_wf == 0
            wf1(:,e)  = temp;
        end
    end

    %% Propagate wavefront
    f_wf1  = fftshift(fft(ifftshift(wf1,1),[],1),1);           % ft wavefront
    I = zeros(length(x),length(E_x),length(z));
%    I_ps = zeros(length(x),length(E_x),length(z));
    for dis = 1:length(z)
%         tic
        f_prop = exp(-1i*pi*lambda.*z(dis).*u.^2); % ft Fresnel propagator
        f_prop = f_prop';
        norm_factor = exp(1i*k.*z(dis))./(2*d_sg1); % for normalation propagator
        norm_factor = norm_factor./abs(norm_factor);
        p_wf = norm_factor.*ifftshift(ifft(f_wf1.*f_prop,[],1),1); % convolute prop with wavefront
        p_I = abs(p_wf).^2; % propagated intensity
        if spherical_wf == 1
            %p_I = p_I;
            % do nothing
        elseif spherical_wf == 0
            % magnification
            M = (d_sg1+z(dis))/d_sg1;
            x_mag = linspace(x(1),x(end),M*N);
            p_I = interp1(x',p_I,x_mag','linear');
            p_I = p_I((round(size(p_I,1)/2)-(N/2)):(round(size(p_I,1)/2)+(N/2)-1),:);
        end
%        I_ps(:,:,dis) = p_I;
        % take source size into account
        w = z(dis)*(source_size/d_sg1); % Consider magnification of source
        prof_source = exp(-(1/2)*(x.^2)/ w^2)'; % gaussian source
        prof_source = prof_source./sum(prof_source); % normalize
        ft_source = fftshift(fft(prof_source,[],1),1); % ft source for convolution
        ft_I = fftshift(fft(p_I,[],1),1); % ft wavefront intensity for convolution
        temp = abs(ifftshift(ifft(ft_I.*ft_source,[],1),1)); % convolution
        I(:,:,dis) = temp;
%         toc
    end
    % % incoherently sum different energy contributions
    I_full = squeeze(sum(I,2));
    clear I
%     I_full_ps = squeeze(sum(I_ps,2)); 
%     clear I_ps

    %% phase stepping
    steps = 15;
    detector_pixsize = 1;
    periods = 2;
    vis = zeros(1,length(z));
    amp = zeros(1,length(z));
    vis1 = zeros(1,length(z));
    amp1 = zeros(1,length(z));
    curves = zeros(steps,length(z));
    %figure, hold on, xlabel('phase step'), ylabel('intensity')
    for dis = 1:length(z)
        M = (d_sg1+z(dis))/d_sg1;
        p2 = p1*M/2;
        [stepping_curve,~,~,~]=phase_stepping(I_full(:,dis),pixsize,detector_pixsize,p2,0.5,steps,periods);
        ft_sc = fft(stepping_curve);
        vis(dis) = abs(ft_sc(1+periods))./abs(ft_sc(1));
        amp(dis) = abs(ft_sc(1));
        vis1(dis) = (max(stepping_curve)-min(stepping_curve))/(max(stepping_curve)+min(stepping_curve));
        amp1(dis) = (max(stepping_curve)+min(stepping_curve))/2;
        curves(:,dis) = stepping_curve;
        %     if mod(pos,10) == 0
    %         plot(stepping_curve, 'o-')
    %     end
    end
    vis_map(:,p) = vis;
    amp_map(:,p) = amp;
    vis_map1(:,p) = vis1;
    amp_map1(:,p) = amp1;
    curves_mat(:,:,p) = curves;
    toc
end

% [PP,ZZ] = meshgrid(p1_range,z); 
% relative_sensitivity = (ZZ.*vis_map.*sqrt(amp_map))./PP;


D_1 = (p1_range.^2)/(8*lambda_from_E(E_0));

d_1 = D_1.*d_sg1./(d_sg1-D_1);
d_1(end) = d_1(end-1);

index = find(p1_range == 7e-06);
d_1(index)

figure, imagesc(vis_map1')
colorbar
yticks(5:5:size(vis_map1,2))
yticklabels(p1_range(5:5:size(vis_map1,2))*1e6)
xticks(50:50:size(vis_map1,1))
xticklabels(z(50:50:size(vis_map1,1)))
xlabel('inter-grating distance [m]')
ylabel('p1 [um]')
title('visibility')
hold on
plot(fliplr(D_1)*length(z)/(z(end)-z(1)),fliplr(1:length(p1_range)),'r-')
plot(fliplr(d_1)*length(z)/(z(end)-z(1)),fliplr(1:length(p1_range)),'r-')



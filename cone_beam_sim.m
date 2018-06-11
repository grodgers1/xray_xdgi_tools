clear all
addpath('/Users/Griffin/Documents/MATLAB/xray_xdgi_tools/tables')
addpath('/Users/Griffin/Documents/MATLAB/xray_xdgi_tools')


p1 = 7e-06; % [m] period of the above mentioned pattern
m1 = 'Au'; % First material (e.g. 'Au')
m2 = 'Au'; % Second material (e.g. 'Au')
t1 = 6e-06; % thickness of first material [m]
t2 = 0; % thickness of second material [m]

source_size = 2e-06; % FWHM of source [m]
d_sg1 = 0.294; % source to g1 distance [m]
z = 0.001:0.001:0.3; % radius of wavefront curvature[m]

x_pixels = 1000; % pixels per period for simulation
reptimes = 20; % repeat g1 so that large z aren't smeared

R = d_sg1;
% Source details
source_spectrum = 0;

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

pixsize = p1/x_pixels;

x = (-((x_pixels*reptimes-1)/2):((x_pixels*reptimes-1)/2))*pixsize;
lambda = lambda_from_E(E_spectrum);
k = 2*pi./lambda;

N = x_pixels*reptimes;
dy0 = reptimes.*(p1/N);
y0 = (-(N/2):(N/2-1)).*dy0;
du = 1./ (N.*dy0); % [1/m] sampling distance in k-space
u  = (-(N/2):(N/2-1)).*du;

g1_pattern = zeros(length(x),length(E_spectrum));
wf1  = zeros(length(x),length(E_spectrum));
for e = 1:length(E_spectrum)
    [delta1,beta1] = get_refindex(m1, E_spectrum(e));
    [delta2,beta2] = get_refindex(m2, E_spectrum(e));
    tmp1 = I_spectrum(e)*exp(-1i*(delta1-1i*beta1)*k(e)*t1);
    tmp2 = I_spectrum(e)*exp(-1i*(delta2-1i*beta2)*k(e)*t2);
    temp = [tmp1*ones(1,round(x_pixels/2)) tmp2*ones(1,round(x_pixels/2))]';
    temp = repmat(temp,reptimes,1);
    g1_pattern(:,e) = temp;
    wf1(:,e)  = temp.*(exp( 1i.*k(e).*sqrt(R.^2+x.^2) ))';
end


% calculates the itensity of a propagated wave field f along the
% z-axis. "psi" is the absolute square root of Eq. (8).
f_wf1  = fftshift(fft(ifftshift(wf1)));           % FFT of wave field
I = zeros(length(x),length(E_spectrum),length(z));
I_ps = zeros(length(x),length(E_spectrum),length(z));
for dis = 1:length(z)
    tic
    
    H   = exp( -1i.*pi.*lambda'.*z(dis).*u.^2); % Fresnel kernel
    H   = H';
    C   = exp(1i.*k'.*z(dis))./(2*R);                    % const
    C   = C./abs(C);                                     % normalized amplitude
    psi = C' .* fftshift(ifft(ifftshift(f_wf1.*H)));   % convolution
    psi = abs(psi).^2;                                % intensity

    % calculates the itensity of a propagated wave field f along the
    % z-axis and taking account of the mutual coherence (principle from
    % Weitkamp-paper with Gaussian src). implements Eqs. (9) and (10).
    w   = z(dis) * (source_size/R);
    sigm_sq  = w^2/(8*log(2)); % FOR FWHM WIDTH?
    srcgauss = exp( -(1/2).* (y0.^2)./ sigm_sq);
    srcgauss = srcgauss./sum(srcgauss);               % normalized Gauss

    gam = fftshift(fft(ifftshift(srcgauss)))';         % damping factor
    Uf  = fftshift(fft(ifftshift(psi)));                % FFT of intensity
    I(:,:,dis)   = abs(fftshift(ifft(ifftshift(Uf.*gam))));
    I_ps(:,:,dis) = psi;
    
    toc
end

I_full = squeeze(sum(I,2));
clear I
I_ps_full = squeeze(sum(I_ps,2)); 
clear I_ps

ptp = 2; % periods to plot
fig_crop = (round(length(x)/2)-(x_pixels*ptp)):(round(length(x)/2)+(x_pixels*ptp));
figure, imagesc(I_full(fig_crop,:)), colormap gray
figure, imagesc(I_full), colormap gray
figure, imagesc(I_ps_full(fig_crop,:)), colormap gray
%figure, imagesc(squeeze(I(fig_crop, round(end/2),:))), colormap gray


p2 = p1;
steps = 15;
detector_pixsize = 1;
periods = 2;
vis = zeros(1,length(z));
%figure, hold on, xlabel('phase step'), ylabel('intensity')
for pos = 1:length(z)
    [stepping_curve,~,~,~]=phase_stepping(I_full(:,pos),pixsize,detector_pixsize,p2,0.5,steps,periods);
    vis(pos) = (max(stepping_curve)-min(stepping_curve))/(max(stepping_curve)+min(stepping_curve));
%     if mod(pos,10) == 0
%         plot(stepping_curve, 'o-')
%     end
end


figure, plot(z,vis)
xlabel('g2 position [m]')
ylabel('visibility')
ylim([0 1])






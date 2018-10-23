clear all

addpath('./functions/')
addpath('./tables/')

%% Naming convention
% l : source-g1
% d : g1-g2
% all distances should be in [m]
% source profile g(x) = 1/(s*sqrt(2*pi)) * exp(-x^2/(2s^2))

%% Setup parameters
E = 50e3; % [eV]
p1 = 7e-06; % [m]
l = 0.4; % [m]
%lpd = 0.296+0.296; % l + d [m]
s = 1.5e-06; % source size [m] 

n = 1; % which talbot order?
nu = 2; % 2 for pi shift, 1 for abs or pi/2

%% Time saving anonymous functions
fd = @(lpd,l) lpd-l;
fM = @(lpd,l) lpd/l;
% parallel beam talbot distance
fDn = @(nu,n,p1,lambda) (1/nu)^2 * (n*p1^2)/(2*lambda);
% cone beam talbot distance
fdn = @(nu,n,p1,lambda,l) l*fDn(nu,n,p1,lambda)/(l-fDn(nu,n,p1,lambda));

%% Functions to find distances
lambda = lambda_from_E(E); % will be in [m]
d = fd(lpd,l);

%% Talbot Distances
Dn = fDn(nu,n,p1,lambda);
dn = fdn(nu,n,p1,lambda,l);

%%
p2 = fM(lpd,l)*p1/nu;

l_range = (30:2.5:50)*1e-2;

dn_l = zeros(size(l_range));
p2_l = zeros(size(l_range));
for i = 1:length(l_range)
    l = l_range(i);
    dn_l(i) = fdn(nu,n,p1,lambda,l);
    
end

figure, plot(l_range, dn_l, 'o')
hold on, plot(l_range,lpd-l_range,'r-')


%% 
% fix lpd = 0.592, fix p2 = 4.8e-6, fix lambda, find p1, l, and dn
E = 50e3;
p2 = 6e-6;
lpd = 0.6; %0.592;

lambda = lambda_from_E(E);

l = lpd/(1+(p2^2/(lpd*2*lambda)))
dn = lpd-l
p1 = 2*l*p2/lpd


%% Visibility reduction from tranverse coherence
% see phd thesis of M. Bech
s = 2e-6;
w = s*dn/l; % demagnified source size [m]
if w/p2 < 1/(2*pi)
    vdt = 1-3.19*(w/p2);
else
    vdt = (8/pi^2)*exp(-2*pi^2*(w/p2)^2);
end


% old settings
s_ak = 2e-6;
w_ak = s_ak*0.296/0.296;
p2_ak = 7e-6;
if w_ak/p2_ak < 1/(2*pi)
    vdt_ak = 1-3.19*(w_ak/p2_ak);
else
    vdt_ak = (8/(pi^2))*exp(-2*(pi^2)*(w_ak/p2_ak)^2);
end



x = 0:0.005:0.5;
x_fwhm = 2.355*x;
figure, plot(x, exp(-(1.887*x_fwhm).^2)) % Timm Weitkamp SPIE
hold on, plot(x,(8/(pi^2))*exp(-2*(pi^2)*x.^2))
vline(w_ak/p2_ak), vline(w/p2)
vline(w_ak/(p2_ak*2), 'b:'), vline(w/(p2*2), 'b:')

%% Calculate grating thickness
% use projection approximation;
hbarc_over_eV = 1.97327e-7; % [m]
trans_g2 = 0.2; % desired transmission of g2

delta = 1.02e-6; % Au @ 55 keV
beta = 1.696e-08; % Au @ 55 keV

t_g1 = pi*hbarc_over_eV/(E*delta);
t_g2 = -log(trans_g2)*hbarc_over_eV/(E*beta);




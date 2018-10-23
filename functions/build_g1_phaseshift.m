function [g1_pattern_phaseshift] = build_g1_phaseshift(m1,m2,t1,t2,dc,E_x,x_pixels,reptimes)
%BUILD_G1 Builds the complex transmission function for the specified
%grating
%   Inputs:
%       m1: material 1 of grating (e.g. 'Au')
%       m2: material 2 of grating (e.g. 'Au')
%       t1: thickness of material 1 [m] (e.g. 6e-6)
%       t2: thickness of material 1 [m] (e.g. 0)
%       dc: duty cycle [0-1]
%               note: if DC = 0.4, this means m2 is wider than m1
%       E_x: x-axis of energy spectrum [eV]
%       E_spectrum: intensity of beam at energies specified by E_x [1]
%       x_pixels: number of pixels per grating period for simulation
%       reptimes: number of repetitions of the grating period
%   Outputs:
%       g1_pattern: complex transmission function of grating, with 2nd
%           dimension corresponding to different energies. For use in
%           simulations using the projection approximation


if dc > 1 || dc < 0
    fprintf('NOTE: Duty cylce must be in [0 1]. Using DC = 0.5 instead. \n')
    dc = 0.5;
end

% total number of pixels
N = x_pixels*reptimes; 
% size of bars
sb1 = ceil(x_pixels*dc); % size of grating bar material 1 (pixels)
sb2 = x_pixels-sb1; % size of grating bar material 2 (pixels)
% useful conversion
lambda = lambda_from_E(E_x)';
k = 2*pi./lambda';
% run loop over energies, build complex transmission function for each
g1_pattern_phaseshift = zeros(N,length(E_x)); % complex transmission function
for e = 1:length(E_x)
    [delta1,~] = get_refindex(m1, E_x(e));
    [delta2,~] = get_refindex(m2, E_x(e));
    tmp1 = delta1*k(e)*t1; % projection approximation
    tmp2 = delta2*k(e)*t2;
    temp = [tmp1*ones(1,sb1) tmp2*ones(1,sb2)]';
    temp = repmat(temp,reptimes,1);
    g1_pattern_phaseshift(:,e) = temp;
end


end


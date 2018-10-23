function [E_spectrum,E_x] = EspectrumGauss(E_0,sig_E,n,E_min,E_max)
%ESPECTRUMGAUSS Makes a gaussian energy spectrum for the source
%   E_spectrum: energy values
%   I_spectrum: source intensity at these energies

E_x = linspace(E_min,E_max,n);
E_spectrum = exp(-(E_x-E_0).^2 / (2*sig_E^2));
% cut energies below 20 keV
E_spectrum(E_x<20e3) = 0;
% normalize
E_spectrum = E_spectrum/sum(E_spectrum);



end


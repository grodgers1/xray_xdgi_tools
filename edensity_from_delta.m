function [rho_e] = edensity_from_delta(delta,E)
%EDENSITY_FROM_DELTA Converts from decrement of Re(index of refraction) to
%electron density. Just to make life easier.
%   delta: decrement of real part of index of refraction (unitless)
%   E: energy [eV]
%   Output is electron density in [electrons per cubic meter]
rho_e = delta.*E^2*1.450492e+27;
end


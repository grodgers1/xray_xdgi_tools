function [ alpha_min ] = min_angle( d,p1,s )
%MIN_ANGLE Summary of this function goes here
%   Detailed explanation goes here
%   All inputs in [m]
%   Take output as relative, not absolute
alpha_min = (p1*0.3./d)*exp((1.877*s*d./(p1*0.3)).^2);
end


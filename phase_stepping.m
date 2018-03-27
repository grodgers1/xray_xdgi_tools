function [ stepping_curve, phase, ac, vis ] = phase_stepping( I, pixsize, detector_pixsize, p2, steps, periods, plot_on )
%phase_stepping This function simulates phase stepping an absorption
% grating over a given interference pattern.
% Inputs:  
%   I_z: interference pattern (intensity)
%   pixsize: pixel for the given intensity pattern [m]
%   detector_pixsize: pixel size of the detector [m] (will implement later)
%   p2: absorption grating period [m]
%   steps: how many phase steps
%   periods: over how many periods do we step?
% Outputs:
%   stepping_curve: the phase stepping curve
%   phase: the phase of this curve (found with fft)
%   ac: absorption value of this curve (found with fft)
%   vis: visibility of this curve (found with fft)

% Make 2nd grating, find stepping positions:
p2_pix = round(p2/pixsize); % size of p2 in pixels
avg_window = 1:round(p2_pix/2); % size of opening in grating, which detector finds avg intensity over
step_positions = 0:round(p2_pix*periods/steps):round(p2_pix*periods); % positions for phase steps
step_positions = step_positions + round(length(I)/2); % do stepping in middle of pattern (avoid blurring from finite height)

% Make stepping curve:
stepping_curve = zeros(1,steps); % initializing the stepping curve
for ps  = 1:steps
    stepping_curve(ps) = mean(I(step_positions(ps)+avg_window));
end

% Get dpc, vis, and absorption signals
ft_stepping_curve = fft(stepping_curve);
ac = abs(ft_stepping_curve(1)); 
vis = abs(ft_stepping_curve(1+periods))./abs(ft_stepping_curve(1));
phase = angle(ft_stepping_curve(1+periods));
if plot_on == 1
   figure, plot(stepping_curve)
end

end


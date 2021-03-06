function [ stepping_curve ] = phase_stepping( I, dx, det_pixsize, p2, duty_cycle, steps, periods )
%phase_stepping This function simulates phase stepping an absorption
% grating over a given interference pattern. It finds the stepping curve
% for one detector pixel of specified size, located at the center of the
% Inputs:  
%   I_z: interference pattern (intensity)
%   dx: pixel size for the given intensity pattern [m]
%   det_pixsize: pixel size of the detector [m]
%   p2: absorption grating period [m]
%   duty_cycle: duty cycle of absorption grating (0-1)
%   steps: numer of phase steps
%   periods: number of periods to step over
% Outputs:
%   stepping_curve: the phase stepping curve

% Make 2nd grating, find stepping positions:
len = length(I);
p2_pix = round(p2/dx); % size of p2 in pixels
det_pix = round(det_pixsize/dx);
if det_pix > len % detector can't be bigger than wavefront
    det_pix = len;
    %fprintf(['Warning: detector pixel bigger than wavefront \n'])
end
reptimes = round(len/p2_pix)+1; % okay to go over, as we cut it to size later

% Make detector
detector = zeros(len,1);
detector(1:det_pix) = ones(det_pix,1);
detector = circshift(detector,round(len/2 - det_pix/2)); % put dector at center

% Make transmission function of second grating
t_p2 = [ones(1,round(p2_pix*duty_cycle)) zeros(1,round(p2_pix*(1-duty_cycle)))];
t_p2 = repmat(t_p2', reptimes,1);

% Set up stepping positions
step_positions = 0:round(p2_pix*periods/steps):round(p2_pix*periods); % positions for phase steps

% Make stepping curve:
stepping_curve = zeros(1,steps); % initializing the stepping curve
for ps  = 1:steps
    temp_g2 = circshift(t_p2, step_positions(ps));
    temp = I.*temp_g2(1:len).*detector;
    stepping_curve(ps) = sum(temp,1);
end

end

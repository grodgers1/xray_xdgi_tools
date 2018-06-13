% p1 in 1e-06 to 10e-6 (1 um to 10 um)
% d in 0.5e-01 to 5.5e-01 (55 cm to 5 cm)

p1 = (1:10)*1e-06; % period of first grating [m]
d = 0.055*(1:10)*1e-01; % Distance from g1 to g2 [m]
s = 4e-06; % source size [m]

alpha_min = zeros(length(p1),length(d))
for p = 1:length(p1)
    for dist = 1:length(d)
        alpha_min(p,dist) = min_angle(d(dist),p1(p),s);
    end
end

figure, imagesc(alpha_min)
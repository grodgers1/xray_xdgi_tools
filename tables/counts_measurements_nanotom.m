% focus 0 has spot size 2 um, focus 1 has spot size 1.5 um

%% Data from counts measurements:
kvp_0 = [10 20 25 30 35 40 45 50 55 60]';
kvp_1 = [30 35 40 45 50 55 60]';

uA_0 = [210 440 550 610 520 350 350 310 340 310]';
uA_1 = [425 380 325 300 275 260 240]';

counts_0 = [0 24 95 210 379 399 559 685 1004 862]';
counts_1 = [90 155 207 272 333 396 455]';

%% Plot
figure, plot(kvp_0, counts_0, 'k*')
hold on, plot(kvp_1, counts_1, 'r*')
xlabel('kV_p')
ylabel('counts')
legend('Focus 0 (2 um)','Focus 1 (1.5 um)')
grid on

%% Estimate mean energy from acceleration voltage
E_mean_0 = ((kvp_0+13)/2);
E_mean_1 = ((kvp_1+13)/2);


figure, plot(kvp_0, E_mean_0)
xlabel('Acceleration voltage (kV_p)')
ylabel('Approx mean energy (keV)')
grid on

%% Fit the counts vs. mean energy
fit_0 = fit(E_mean_0,counts_0,'poly2');
fit_1 = fit(E_mean_1,counts_1,'poly2');

figure, plot(fit_0,E_mean_0,counts_0, 'k*')
hold on, plot(fit_1,E_mean_1,counts_1, 'r*')

counts_focus0 = @(E_mean) (fit_0.p1*E_mean.^2 + fit_0.p2*E_mean + fit_0.p3);
counts_focus1 = @(E_mean) (fit_1.p1*E_mean.^2 + fit_1.p2*E_mean + fit_1.p3);




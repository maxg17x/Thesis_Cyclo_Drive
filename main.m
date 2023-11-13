clear; clc; clf;
geo_param = [33.5,4.5,15,1.45, 0, 0.025];
mat_param = [0.28,0.28,190*10^9,190*10^9,2000,2000,10,0.002];
X = CDA(geo_param, mat_param);
%plot_eff_vs_torque_wOff();
plot_eff_vs_torque_wEcc();
%plot_force_vs_eccentricity();
%plot_backlash_vs_eccentricity();
%plot_eff_vs_ratio()
%X.efficiency(50);
%X.F_cs\
%X.theta_shift(pi/2);
%plot_force_vs_offset();
%X.plot_stiffness();
% X.hertz_stiffness();
% X.pin_backlash();
tic;
% X.pin_backlash();
% X.hertz_stiffness();
% k_vals = X.K
% X.force_distribution();





function plot_eff_vs_torque_wEcc(X)
    eccs = [0.75 1 1.25 1.5 1.75]; test_loads = [1:1:200];
    mat_param = [0.28,0.28,190*10^9,190*10^9,2000,2000,10,0.002];
    eff_array = [];
    hp = {};
    figure(9);
    for i = 1:length(eccs)
        arr = zeros(1,length(test_loads));
        geo_param = [33.5,4.5,15,eccs(i), 0, 0.02];
        X = CDA(geo_param, mat_param);
        X.theta_shift(pi/2);
        for i = 1:length(test_loads)
            arr(i) = X.efficiency(test_loads(i));
        end
        eff_array = [eff_array; arr];
        figure(9);
        hp{i} = plot(1:200,arr); hold on;
    end
    xlim([0 200]);
    xlabel('Output Torque (Nm)'); ylabel('Efficiency');
    hleg = legend('0.75 mm','1.00 mm', '1.25 mm', '1.50 mm', '1.75 mm');
    set(get(hleg,'Title'),'String','Eccentricity (mm)')
end
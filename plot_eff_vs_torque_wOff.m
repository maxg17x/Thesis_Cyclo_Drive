function plot_eff_vs_torque_wOff(X)
    offs = [0.005 0.01 0.02 0.05]; test_loads = [1:1:200];
    mat_param = [0.28,0.28,190*10^9,190*10^9,2000,2000,10,0.002];
    off_array = [];
    hp = {};
    figure(9);
    for i = 1:length(offs)
        arr = zeros(1,length(test_loads));
        geo_param = [33.5,4.5,15,1.2, 0, offs(i)];
        X = CDA(geo_param, mat_param);
        X.theta_shift(pi/2);
        for i = 1:length(test_loads)
            arr(i) = X.efficiency(test_loads(i));
        end
        [~, idx] = min(arr);
        arr(1:idx) = NaN;
        off_array = [off_array; arr];
        figure(9);
        hp{i} = plot(1:200,arr); hold on;
    end
    xlim([0 200]);
    xlabel('Output Torque (Nm)'); ylabel('Efficiency');
    hleg = legend('0.005 mm','0.010 mm', '0.020 mm', '0.050 mm');
    set(get(hleg,'Title'),'String','Offset (mm)')
end
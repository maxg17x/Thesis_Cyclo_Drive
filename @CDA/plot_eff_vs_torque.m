function plot_eff_vs_torque(X)
    test_loads = [0:1:200];
    eff_array = zeros(1,length(test_loads));
    eccs = [0.75 1 1.25 1.5 1.75];
    figure(6);
    for i = 1:length(test_loads)
        eff_array(i) = X.efficiency(test_loads(i));
    end
    plot(0:1:200,eff_array(1,:)); hold on;
    xlim([0 200]); ylim([0.8 1]); grid on;
    xlabel('Output Torque (Nm)'); ylabel('Efficiency');
end
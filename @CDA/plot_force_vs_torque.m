function plot_force_vs_torque(X)
    test_loads = [10, 25, 50, 100, 150, 200];
    force_array = [];
    figure(6);
    for i = 1:length(test_loads)
        X.force_distribution(test_loads(i));
        force_array = [force_array; X.F_c];
        %force_array(force_array == 0) = NaN;
        hp{i} = plot(1:(X.N+1),0.001*force_array(i,:)); hold on;
        plot(1:(X.N+1),0.001*force_array(i,:),'k.'); hold on;
    end
    xlim([3 9]); ylim([0 4.5]); grid on;
    xlabel('Pin Number'); ylabel('Force on pin (kN)');
    hleg = legend([hp{1}(1);hp{2}(1);hp{3}(1);hp{4}(1);hp{5}(1);hp{6}(1)], ...
    '10 Nm','25 Nm', '50 Nm', '100 Nm', '200 Nm', '400 Nm');
    set(get(hleg,'Title'),'String','Output Torque (Nm)')
end
function plot_force_vs_offset(X)
    offsets = [0.005 0.01 0.02 0.05 0.1];
    mat_param = [0.28,0.28,190*10^9,190*10^9,2000,2000,10,0.002];
    force_array = [];
    figure(7);
    for i = 1:length(offsets)
        geo_param = [33.5,4.5,15,1.45, 0, offsets(i)];
        X = CDA(geo_param, mat_param);
        X.theta_shift(pi/2);
        X.force_distribution(100);
        force_array = [force_array; X.F_c];
        %force_array(force_array == 0) = NaN;
        figure(7);
        hp{i} = plot(1:(X.N+1),0.001*force_array(i,:)); hold on;
        plot(1:(X.N+1),0.001*force_array(i,:),'k.'); hold on;
    end
    xlim([3 9]); ylim([0 4.5]); grid on;
    xlabel('Pin Number'); ylabel('Force on pin (kN)');
    hleg = legend([hp{1}(1);hp{2}(1);hp{3}(1);hp{4}(1);hp{5}(1)], ...
    '0.005 mm','0.01 mm', '0.02 mm', '0.05 mm', '0.1mm');
    set(get(hleg,'Title'),'String','Offset Modification (mm)')
end
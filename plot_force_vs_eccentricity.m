function plot_force_vs_eccentricity(X)
    eccs = [0.75 1 1.25 1.5 1.75];
    mat_param = [0.28,0.28,190*10^9,190*10^9,2000,2000,10,0.002];
    force_array = [];
    figure(8);
    for i = 1:length(eccs)
        geo_param = [33.5,4.5,15,eccs(i), 0, 0.02];
        X = CDA(geo_param, mat_param);
        X.theta_shift(pi/2);
        X.force_distribution(100);
        force_array = [force_array; X.F_c];
        %force_array(force_array == 0) = NaN;
        figure(8);
        hp{i} = plot(1:(X.N+1),0.001*force_array(i,:)); hold on;
        plot(1:(X.N+1),0.001*force_array(i,:),'k.'); hold on;
    end
    xlim([3 9]); ylim([0 7]); grid on;
    xlabel('Pin Number'); ylabel('Force on pin (kN)');
    hleg = legend([hp{1}(1);hp{2}(1);hp{3}(1);hp{4}(1);hp{5}(1)], ...
    '0.75 mm','1.00 mm', '1.25 mm', '1.50 mm', '1.75 mm');
    set(get(hleg,'Title'),'String','Eccentricity (mm)')
end
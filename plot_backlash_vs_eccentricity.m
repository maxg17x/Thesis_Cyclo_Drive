function plot_backlash_vs_eccentricity()
    offs = [0.005 0.01 0.02 0.05 0.1];
    mat_param = [0.28,0.28,190*10^9,190*10^9,2000,2000,10,0.002];
    for j = 1:5
        geo_param = [36,4.5,15,1.45, 0, offs(j)];
        X = CDA(geo_param, mat_param);
        res = 360*1.5;
        bkl_angles = zeros(X.N+1, res);
        for i = 1:res
            theta_shift(X, i*3*pi/res);
            bkl_angles(1:(X.N+1), i) = X.bkl_angle;
        end
        bkl_angles(isinf(bkl_angles)) = NaN
        figure(4);
        hp{j} = plot(1:res,rad2deg(min(bkl_angles(:,:)))); hold on;
        xlim([45 275]); ylim([0.03 0.12]);
        
    end
    xlabel('Input Angle (deg)'); ylabel('Backlash (deg)');
    hleg = legend([hp{1}(1);hp{2}(1);hp{3}(1);hp{4}(1);hp{5}(1)],'0.005 mm','0.01 mm', '0.02 mm', '0.05 mm', '0.1 mm');
    set(get(hleg,'Title'),'String','Offset Mod. (mm)')
end
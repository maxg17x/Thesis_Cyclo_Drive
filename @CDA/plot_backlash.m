function plot_backlash(X)
    figure(2);
    res = 360*1.5;
    bkl_angles = zeros(X.N+1, res);
    for i = 1:res
        theta_shift(X, i*3*pi/res);
        bkl_angles(1:(X.N+1), i) = X.bkl_angle;
    end
    bkl_angles(isinf(bkl_angles)) = NaN
    for i = 1:(X.N+1)
        plot(1:res,rad2deg(bkl_angles(i,:))); hold on;
    end
    plot(1:res,rad2deg(min(bkl_angles(:,:))),'r','LineWidth',3); hold on;
    xlim([0 320]); ylim([0.05 0.06]);
    xlabel('Input Angle (deg)'); ylabel('Static Backlash (deg)');
end
function plot_stiffness(X)
    figure(3);
    res = 360*1.5;
    K_values = zeros(X.N+1, res);
    for i = 1:res
        theta_shift(X, i*3*pi/res);
        X.hertz_stiffness();
        K_values(1:(X.N+1), i) = X.K;
    end
    for i = 1:(X.N+1)
        figure(3);
        plot(1:res,K_values(1,:)); hold on;
    end
    xlim([0 360]);
     xlabel('Input Angle (deg)'); ylabel('Stiffness Constant (N/m)');
end
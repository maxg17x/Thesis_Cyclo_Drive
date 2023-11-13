function theta_shift(X, t_new)
    X.t = t_new; X.ecc_pos = X.E*[cos(X.t); sin(X.t)];
    X.profile_t = frame_trans(X, X.profile_b, X.ecc_pos, -X.t/X.N);
    set(X.h_profile,'xdata',X.profile_t(1,:),'ydata',X.profile_t(2,:));
    set(X.h_center_pt, 'xdata', X.E*cos(X.t), 'ydata', X.E*sin(X.t));
    X.pin_backlash();
    for i = 1:length(X.cont_angle)
        if X.cont_angle(i) > 2*pi/X.N
            X.cont_angle(i) = X.cont_angle(i)-(2*pi/X.N);
        end
    end

    refreshdata
    drawnow
end
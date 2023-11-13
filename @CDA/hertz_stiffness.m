% Calculate the approximated hertzian stiffness constant for each point
function hertz_stiffness(X)
    X.K = zeros(1,length(X.bkl_angle));
    for i = 1:(X.N+1)
        if isinf(X.bkl_arc(i))
            continue
        end
        % Find r_ci, the approxinmate radius of the gear profile at cont_pt
        fit_span = 0.0005;
        p1 = base_eqn(X,X.dr_e,X.dr_p,X.cont_angle(i)-fit_span);
        p2 = base_eqn(X,X.dr_e,X.dr_p,X.cont_angle(i));
        p3 = base_eqn(X,X.dr_e,X.dr_p,X.cont_angle(i)+fit_span);
        fit_pts = frame_trans(X,[p1 p2 p3],X.ecc_pos, -X.t/X.N);
        % plot(fit_pts(1,:),fit_pts(2,:),'k.','MarkerSize',8); hold on;
        [r_ci, xy_i] = fit_circle_through_3_points(fit_pts.');
        ang = 0:0.01:2*pi;
        % plot(xy_i(1)+r_ci*cos(ang),xy_i(2)+r_ci*sin(ang))
        % Determine relative concavity of the fitted circle - CHECK ME
        mid_r = (X.r_e + X.dr_e) - (X.r_p + X.dr_p);
        % plot(X.ecc_pos(1)+mid_r*cos(ang),X.ecc_pos(2)+mid_r*sin(ang))
        if norm(X.cont_pt(i)-X.ecc_pos) < (X.r_e-X.r_p)
            r_ci = -r_ci;
        end
        % 
        % Linearise the stiffness hertzian stiffness relationship
            % Calculate the equivalent radius
        rho_star = (10^-3)*norm(r_ci*X.r_p/(r_ci+X.r_p));
            % Calculate the contact-rectangle half-width for x_fit
        x_fit = X.LR_p/2;
        b_fit = sqrt(4*x_fit*rho_star/(pi*X.B*0.001*X.E_star));
        exp1 = ((1-X.v_c^2)/X.E_c)*(log(4*norm(r_ci)*0.001/b_fit)-0.5);
        exp2 = ((1-X.v_p^2)/X.E_p)*(log(4*X.r_p*0.001/b_fit)-0.5);
        X.K(i) = 0.5*pi*X.B*0.001/(exp1 + exp2); % N/m
    end
end
% Calculate the force distribution using bkl_arc length and stiffness values
function force_distribution(X, max_torque)
    % Determine the normal vectors for contact
    d_c = zeros(2,X.N+1);
    r_c = zeros(2,X.N+1);
    contact_lines = [];
    moment_lines = [];
    for i = 1:(X.N+1)
        contact_vector = (X.pin_pos(:,i) - X.cont_pin_pt(:,i));
        d_c(:,i) = contact_vector/norm(contact_vector);
        r_c(:,i) = X.cont_pin_pt(:,i) - X.ecc_pos;
        contact_lines = [contact_lines,X.pin_pos(:,i),X.pin_pos(:,i)-X.r_p*d_c(:,i),[nan; nan]];
        moment_lines = [moment_lines,X.ecc_pos,X.cont_pin_pt(:,i),[nan; nan]];
    end
    set(X.h_cont_lines,'xdata',contact_lines(1,:),'ydata',contact_lines(2,:));
    set(X.h_moment_lines,'xdata',moment_lines(1,:),'ydata',moment_lines(2,:));
    dir_const = 0.001*(d_c(1,:) .* r_c(2,:) - d_c(2,:) .* r_c(1,:));
    dir_const(dir_const > 0 | isnan(dir_const)) = 0;
    % Sort pins in order of backlash angle, retain the indexes
    finiteIndices = find(isfinite(X.bkl_angle));
    [~, sortIndices] = sort(X.bkl_angle(finiteIndices));
    pin_contact_order = finiteIndices(sortIndices);
    % Determine the number of pins needed to resist the torque and their force
    for i = 1:length(pin_contact_order)
        syms d;
        pin_torques = sym(zeros(1,i));
        for j = 1:i
            idx = pin_contact_order(j);
            pin_torques(j) = dir_const(idx)*X.K(idx)*(-d + 0.001*X.bkl_arc(idx));
        end
        y = vpa(pin_torques, 4);
        dist_needed = double(solve(max_torque == sum(pin_torques), d));
        bkl_max = 0.001*X.bkl_arc(pin_contact_order(i+1));
        if dist_needed > 0.001*X.bkl_arc(pin_contact_order(i+1))
            continue
        else
            idc = pin_contact_order;
            D = dir_const;
            K = X.K;
            x = 0.001*X.bkl_arc;
            gamma = dist_needed;
            X.F_c = X.K.*(gamma - 0.001*X.bkl_arc);
            X.F_c(X.F_c < 0) = 0;
            break
        end
    end
end
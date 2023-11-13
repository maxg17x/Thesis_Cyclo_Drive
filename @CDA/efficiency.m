function eff = efficiency(X, max_torque)
    X.force_distribution(max_torque);
    % Get the force of friction at the contact points
    F_friction = X.u*X.F_c;
    % Get the relative velocity vector at each contact point
    point_velocity = zeros(3,X.N+1);
    for i = 1:(X.N+1)
        % Net velocity at contact pt
        point_velocity(:,i) = cross([0; 0; X.w*pi/30], [0.001*X.ecc_pos; 0]) + ...
            cross([0; 0; -(1/(X.N+1))*X.w*pi/30], 0.001*[X.cont_pin_pt(:,i) - X.ecc_pos; 0]);
        % Velocity in the direction of friction
        contact_vector = (X.pin_pos(:,i) - X.cont_pin_pt(:,i));
        contact_vector = contact_vector/norm(contact_vector);
        n = [-contact_vector(2); contact_vector(1); 0];
        point_velocity(:,i) = (dot(point_velocity(:,i),n) / dot(n, n)) * n;
    end
    P_losses = F_friction.*vecnorm(point_velocity);
    P_losses(isnan(P_losses)) = 0;
    s = sum(P_losses)
    P_out = max_torque*(X.w/X.N)*pi/30;
    eff = (P_out-(max_torque*(X.N-1)/(X.N)))/(P_out+s);
end
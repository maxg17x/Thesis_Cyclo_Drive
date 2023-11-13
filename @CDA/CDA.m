classdef CDA < handle
properties
    % Geometric Parameters
    r_e         % Pitch circle radius of external pins
    r_p         % Radius of external pins
    N           % Reduction ratio
    E           % Shaft eccentricity
    dr_e        % Isometric modifier (delta r_e)
    dr_p        % Offset modifier (delta r_p)
    % Material Parameters
    v_c         % Poisson's ratio for cycloidal gear
    v_p         % Poisson's ratio for external pin
    E_c         % Elastic modulus for cycloidal gear
    E_p         % Elastic modulus for external pin
    E_star      % Equivalent hertzian elastic modulus
    B           % Width of cycloidal gear
    LR_c        % Static load rating for cycloidal gear
    LR_p        % Static load rating for external pin
    u           % Coefficient of friction for external pins / cyclo
    % Dynamic Parameters
    w           % Shaft input speed
    t           % Shaft input angle
    profile_b   % Base cycloidal gear equation
    profile_t   % Cycloidal gear equation for current t
    % Construction Parameters
    pin_pos     % Array of external pin positions
    ecc_pos     % Position of eccentric centre given t
    % Backlash Variables
    cont_angle  % Projected contact point angle array
    cont_pt     % Projected contatct points on cyclodial gear
    cont_pin_pt % Projected contact points on the external pins
    bkl_arc     % Backlash arc length for external pins
    bkl_angle   % Backlash angle array for external pins    
    % Force Variables
    d_c         % External pin contact unit vectors (for direction)
    K           % Hertzian stiffness values for each contacting pin
    F_c         % External pin contact forces
    % Handles
    h_profile
    h_cont_pt
    h_center_pt
    h_cont_lines
    h_moment_lines
    
end
methods

    % Assign varaibles, setup figures, plot the graphical representation
    function X = CDA(geo_param, mat_param)
        % Extract geometric and material parameters
        geo_cell = num2cell(geo_param); mat_cell = num2cell(mat_param);
        [X.r_e,X.r_p,X.N,X.E,X.dr_e,X.dr_p]         = geo_cell{:};
        [X.v_c,X.v_p,X.E_c,X.E_p,X.LR_c,X.LR_p,X.B,X.u] = mat_cell{:};
        X.t = 0; X.w = 1000; X.ecc_pos = X.E*[cos(X.t); sin(X.t)];
        % Define the poin position as an array of coordinates {2 by n}
        X.pin_pos = (X.r_e)*[cos(2*pi*(1:(X.N+1))/(X.N+1));
            sin(2*pi*(1:(X.N+1))/(X.N+1))];
        % Create the handles for animated graphs
        figure(1);
        X.h_profile = plot(NaN,NaN,'m','LineWidth',1); hold on;
        X.h_cont_pt = plot(0,0,'r.','MarkerSize',6); hold on;
        X.h_center_pt = plot(X.E*cos(X.t), X.E*sin(X.t), 'm.', 'Markersize', 6); hold on;
        X.h_cont_lines = plot(NaN,NaN,'g','LineWidth',2); hold on;
        X.h_moment_lines = plot(NaN,NaN,'--b','LineWidth',0.2); hold on;
        % Generate the cycloidal gear base profile with initial conditions
        profile_res = 5000; n = linspace(0,2*pi,profile_res);
        X.profile_b = frame_trans(X,base_eqn(X,X.dr_e,X.dr_p,n),0,X.t);
        X.profile_t = frame_trans(X,base_eqn(X,X.dr_e,X.dr_p,n), X.ecc_pos,X.t);
        set(X.h_profile,'xdata',X.profile_t(1,:),'ydata',X.profile_b(2,:));
        % Generate the exterior pin graphs, global centre
        k = linspace(0,2*pi,1000); figure(1);
        for i = 1:(X.N+1)
            xi = X.pin_pos(1,i)+X.r_p*cos(k); yi = X.pin_pos(2,i)+X.r_p*sin(k);
            plot(xi,yi,'k','LineWidth',1); hold on; axis equal;
        end
        plot(0,0,'k.','MarkerSize',8); hold on;
        xlim([-45 45]); ylim([-45 45]);
        % Calculate the hertzian reduced elastic modulus
        X.E_star = X.E_c*X.E_p/((1-X.v_c^2)*X.E_p+(1-X.v_p^2)*X.E_c);
        X.pin_backlash();
        X.hertz_stiffness();
        % X.force_distribution();
    end

    % Define the base equation of the cycloidal profile
    x = base_eqn(X, dr_e, dr_p, t);
    
    % Homogenous transformation function
    out = frame_trans(~, vec, shift, angle);

    % Modify the graph and parameters for a new t (angle) value
    theta_shift(X, t_new);

    % For the current instant t (angle) calcaulate the pin backlash angles
    pin_backlash(X);

    % Find the projected backlash for given contact angle estimate   
    angle = proj_angle(X,i,c_angle);

    % Find the closest contact point given the current estimate
    contact_finder(X, i, res);

    % Calculate the approximated hertzian stiffness constant for each point
    hertz_stiffness(X)  

    % Helper function for fitting circles for hertxian contact
    [R,xcyc] = fit_circle_through_3_points(~, ABC)

    % Calculate the force distribution using bkl_arc length and stiffness values
    force_distribution(X, max_torque);

    % Calculate the efficiency of the drive given some max torque
    eff = efficiency(X, max_torque);

    % Plot the equivalent backlash for all pins
    plot_backlash(X);

    % Plot the stiffness constant vs. angle graph for one pin
    plot_stiffness(X);

    % Display a graph of the contact forces on each pin
    plot_force_distribution(X);

    plot_eff_vs_torque(X);
end
end
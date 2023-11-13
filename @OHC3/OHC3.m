classdef OHC3 < handle
properties
    % Input parameters
    r_e; r_p; N; E; w;                                                    % PCR of ext. pins, pin radius, reduction ratio, shaft eccentricity
    dr_e; dr_p;                                                        % Isometric modification, Offset modification
    theta; E_red;                                                      % Current input angle, reduced elastic modulus
    % HC profiles
    profile_base_fixed; profile_base_free;                             % HC profile without modifications at theta = 0 and current theta
    profile_fixed; profile_free;                                       % HC profile with modifcations at theta = 0 and current theta
    % Other properties
    pin_i_base; pin_i_offset; p_e                                      % Point vectors for the locations of pin k, eccentric center
    first_cont_pin_num; p_c; Pp_c; cont_angles;                              % First contacting pin, array of contact location vectors, contact angles
    est_cont_angle; angle_cont;
    Fc; hat_Fc; gamma_g                                                % Store an array of contact force vectors, direction vectors and backlash
    h_current; h_p_e; h_cont_pt; h_backlash; h_force_distribution;         % Handles for dynamically plotting
    h_no_contact; h_zero_line_c, h_zero_line_nc; h_speed_line;
    % Buffers
    backl_theta;
end
methods
    function X = OHC3(r_e, r_p, E, N, dr_e, dr_p, E_c, v_c, E_p, v_p, w)   % Initialisation function
        X.r_e = r_e; X.r_p = r_p; X.E = E; X.N = N; X.dr_e = dr_e; X.w = w;  % Set all given variables
        X.dr_p = dr_p; X.theta = 0; X.p_e = E*[cos(X.theta); sin(X.theta)];
        X.E_red = E_c*E_p/((1-v_c^2)*E_p+(1-v_p^2)*E_c);               % Calculate the reduced elastic modulus                                               % Generate initial profiles for theta = 0 and pin vectors;        
        X.backl_theta = []; X.Pp_c = zeros(2,N+1);
        X.profile_base_fixed = generate_profile(X,0,0,0); 
        X.profile_base_free  = generate_profile(X,0,0,X.p_e);
        X.profile_fixed      = generate_profile(X,dr_e,dr_p,0);
        X.profile_free       = generate_profile(X,dr_e,dr_p,X.p_e);
        X.pin_i_base   = r_e*[cos(2*pi*(1:(N+1))/(N+1)); sin(2*pi*(1:(N+1))/(N+1))];
        X.pin_i_offset = (r_e+dr_e)*[cos(2*pi*(1:(N+1))/(N+1)); sin(2*pi*(1:(N+1))/(N+1))];
        % Plot the exterior pins, centre and profiles
        k = linspace(0,2*pi,100);
        figure(1);
        for i = 1:(X.N+1)
            xi = X.pin_i_base(1,i)+r_p*cos(k); yi = X.pin_i_base(2,i)+r_p*sin(k);
            plot(xi,yi,'k','LineWidth',1); hold on; axis equal;
        end
        plot(0,0,'k.','MarkerSize',8); hold on;
        X.h_current = plot(0,0,'m','LineWidth',1); hold on;
        X.h_zero_line_c = plot(0,0,'m--','LineWidth',0.25); hold on;
        X.h_no_contact = plot(0,0,'b','LineWidth',0.25); hold on;
        X.h_zero_line_nc = plot(0,0,'b--','LineWidth',0.25); hold on;
        X.h_p_e = plot(X.p_e(1),X.p_e(1),'b.','MarkerSize',8); hold on;
        X.h_cont_pt = plot(0,0,'r.','MarkerSize',6); hold on;
        X.h_backlash = plot(0,0,'m','LineWidth',1); hold on;
        X.h_speed_line = plot(0,0,'g','LineWidth',1); hold on;
        X.h_force_distribution = plot(0,0,'r','LineWidth',1.5); hold on;
    end

    % ----------------------------------------------------------------------------------------------------------------------------------

    function profile = generate_profile(X,dr_e,dr_p,shift)
        profile_res = 5000; t = linspace(0,2*pi,profile_res);
        profile = trans(X,base_eqn(X,dr_e,dr_p,t),shift,-X.theta/X.N);     % Apply a linear transformation to get the profile at current theta
    end
        
    % ----------------------------------------------------------------------------------------------------------------------------------
        
    function theta_shift(X, theta)
        X.theta = theta; X.p_e = X.E*[cos(theta); sin(theta)];             % Update theta and the eccentric center point vector
        X.profile_base_free = trans(X,X.profile_base_fixed,X.p_e,-X.theta/X.N); % Use the fixed profiles to obtain the profiles for current theta
        X.profile_free      = trans(X,X.profile_fixed,X.p_e,-X.theta/X.N);
        set(X.h_no_contact,'xdata',X.profile_free(1,:),'ydata',X.profile_free(2,:));
        set(X.h_p_e,'xdata',X.p_e(1),'ydata',X.p_e(2));
        X.pin_backlash();
        % Normalise bad backlash angles
        for i = 1:length(X.angle_cont)
            if X.angle_cont(i) > 2*pi/X.N
                X.angle_cont(i) = X.angle_cont(i)-(2*pi/X.N);
            end
        end
        %
        X.profile_free = trans(X,X.profile_fixed,X.p_e,min(X.angle_cont)-X.theta/X.N);
        set(X.h_current,'xdata',X.profile_free(1,:),'ydata',X.profile_free(2,:));
        X.backl_theta = [X.backl_theta [X.theta; X.angle_cont.']];
        contact_lines = [];
        for i = 1:(X.N+1)
            if isinf(X.angle_cont(i))
                X.hat_Fc(:,i) = [0; 0]; 
                cont_vec = [0; 0];
            else
                cont_vec = X.pin_i_base(:,i)-X.Pp_c(:,i);
                X.hat_Fc(:,i) = cont_vec/norm(cont_vec);
                
            end
            contact_lines = [contact_lines,X.pin_i_base(:,i),X.pin_i_base(:,i)-X.r_p*cont_vec/norm(cont_vec),[nan; nan]];
        end
        set(X.h_force_distribution,'xdata',contact_lines(1,:),'ydata',contact_lines(2,:));
        % Produce the plot of the zero lines for nc and c OHCs
        pt_zero_nc = trans(X,base_eqn(X,X.dr_e,X.dr_p,0),X.p_e,-X.theta/X.N);
        set(X.h_zero_line_nc,'xdata',[X.p_e(1); pt_zero_nc(1)],'ydata',[X.p_e(2); pt_zero_nc(2)]); 
        pt_zero_c = trans(X,base_eqn(X,X.dr_e,X.dr_p,0),X.p_e,min(X.angle_cont)-X.theta/X.N);
        set(X.h_zero_line_c,'xdata',[X.p_e(1); pt_zero_c(1)],'ydata',[X.p_e(2); pt_zero_c(2)]);
        % Hertz stiffness cacls
%         k_array = zeros(1,X.N+1);
%         for i = 1:(X.N+1)
%             if isinf(X.angle_cont(i))
%                 k_array(i) = 0;
%                 continue
%             else
%                 k_array(i) = hertz_stiffness(X,i);
%             end
%         end
%         k_array;
          % SPEED CALCS
%         ptx  = contact_lines(:,8);
%         v_c2 = contact_speed(X,ptx);
%         speed = norm(v_c2)
%         pt1 = ptx; pt2 = pt1 + v_c2/5;
%         set(X.h_speed_line,'xdata',[pt1(1); pt2(1)],'ydata',[pt1(2); pt2(2)])
        refreshdata
        drawnow
    end

    % ----------------------------------------------------------------------------------------------------------------------------------

    function pin_backlash(X)                                               % Iteratively calculate backlash using smallest arc lenth
        X.est_cont_angle = pi/X.N + (0:1:X.N).*(2*pi/X.N);                 % Initial guesses for closest contact points are at the peaks of OHC
        X.p_c = trans(X,base_eqn(X,X.dr_e,X.dr_p,X.est_cont_angle),X.p_e,-X.theta/X.N);
        start_res = 0.01;                                                  % Angular sampling resolution
        for i = 1:(X.N+1)
            X.angle_cont(i) = X.proj_angle(i,X.est_cont_angle(i));
        end
        X.angle_cont(isnan(X.angle_cont)) = Inf;
        for i = 1:length(X.est_cont_angle)                                     
            if isinf(X.est_cont_angle(i))                                  % Skip calculations for pins that will never be contacted
                continue
            end
            X.contact_finder(i,start_res)                               % Recursively calculate the contact point and angle
        end
    end

    % ----------------------------------------------------------------------------------------------------------------------------------

    function angle = proj_angle(X,i,c_angle)                               % Find the projected backlash for given contact angle estimate                    
        temp_p_ci = trans(X,base_eqn(X,X.dr_e,X.dr_p,c_angle),X.p_e,-X.theta/X.N);                                             
        r1 = norm(temp_p_ci - X.p_e); c1 = X.p_e;                          % Find the point on the OHC given some c_angle
        r2 = X.r_p; c2 = X.pin_i_base(:,i);                                % project an arc from the contact angle estimate to the pin
        [x,y] = circcirc(c1(1),c1(2),r1,c2(1),c2(2),r2); p = [0; 0];
        if norm(temp_p_ci-[x(1);y(1)]) <= norm(temp_p_ci-[x(2);y(2)])      % Find the intersection point between the pin and the arc
            p = [x(1);y(1)];
        else
            p = [x(2);y(2)];
        end                                                                
        X.Pp_c(:,i) = p;
        %plot(p(1),p(2),'g*','MarkerSize',6);
        v1 = temp_p_ci-X.p_e; v2 = p-X.p_e;                                % Vectors that subtend the arc of intrest
        angle = acos(dot(v1,v2)/(norm(v1)*norm(v2)));                      % Calculate the angle of the projected arc
        angle(isnan(angle)) = Inf;                                         % nan -> Inf filter                                        
    end

    % ----------------------------------------------------------------------------------------------------------------------------------

    function contact_finder(X,i,res)                                       % Find the closes contact point given the current estimates
        X.p_c(:,i) = trans(X,base_eqn(X,X.dr_e,X.dr_p,X.est_cont_angle(i)),X.p_e,-X.theta/X.N);
        set(X.h_cont_pt,'xdata',X.p_c(1,:),'ydata',X.p_c(2,:));
        % plot(X.p_c(1,i),X.p_c(2,i),'m*','MarkerSize',6);
        res_thresh = 10^(-5);                                              % Recursively improve estimates till res_thresh is reached
        if res > res_thresh
            temp = X.est_cont_angle(i);                                    % Check if the left or right points offer improvements
            left_in = temp+res; right_in = temp-res; curr_out = X.angle_cont(i);
            left_out = X.proj_angle(i,left_in); right_out = X.proj_angle(i,right_in);
            if left_out < curr_out
                X.est_cont_angle(i) = left_in;
                X.angle_cont(i) = left_out;
                X.contact_finder(i,res);
            elseif right_out < curr_out
                X.est_cont_angle(i) = right_in;
                X.angle_cont(i) = right_out;
                X.contact_finder(i,res);
            else
                X.contact_finder(i,0.5*res);                            % If the left/right estimates are worse lower search resolution
            end
        end
    end

    % ----------------------------------------------------------------------------------------------------------------------------------

    function plot_backlash(X)
        figure(2);
        plot(rad2deg(X.backl_theta(1,:)),rad2deg(X.backl_theta(2:end,:))); hold on;
        plot(rad2deg(X.backl_theta(1,:)),rad2deg(min(X.backl_theta(2:end,:))),'LineWidth',3);
        xlim([0 360]); ylim([0 4]);
        xlabel('Input Angle (deg)'); ylabel('Static Backlash');
    end

    % ----------------------------------------------------------------------------------------------------------------------------------

    function ki = hertz_stiffness(X,i)                                  % Calculate the approximated hertzian stiffness constant
        % Find r_ci, the approximate radius of the OHC at the contact point
        fit_span = 0.001; syms x y; fit_pts = zeros(2,3);                  % Set the gap between fitting points
        fit_pts(:,1) = trans(X,base_eqn(X,X.dr_e,X.dr_p,X.est_cont_angle(i)-fit_span),X.p_e,-X.theta/X.N);
        fit_pts(:,2) = trans(X,base_eqn(X,X.dr_e,X.dr_p,X.est_cont_angle(i)),X.p_e,-X.theta/X.N);
        fit_pts(:,3) = trans(X,base_eqn(X,X.dr_e,X.dr_p,X.est_cont_angle(i)+fit_span),X.p_e,-X.theta/X.N);% Points for the OHC contactt arc to be fitted with
        % Plot the fit points for verification
        %plot(fit_pts(1,:),fit_pts(2,:),'k.','MarkerSize',8); hold on;
        syms x y r;
        eq1 = (fit_pts(1,1)-x)^2+(fit_pts(2,1)-y)^2-r^2==0;
        eq2 = (fit_pts(1,2)-x)^2+(fit_pts(2,2)-y)^2-r^2==0;
        eq3 = (fit_pts(1,3)-x)^2+(fit_pts(2,3)-y)^2-r^2==0;
        [x_c, y_c] = solve([eq2-eq1,eq3-eq1],[x,y]);
        r_ci = solve(subs(eq1,[x,y],[x_c,y_c]),r);                      % Solve the system of equations to obtain r_ci
        a = X.p_c(:,i)-X.p_e; b = [x_c; y_c]-X.p_c(:,i);
        % Plot the circle for verification
        %t = linspace(0,2*pi,500); c = r_ci(2)*[cos(t); sin(t)];
        %plot(x_c+c(1,:),y_c+c(2,:),'k','LineWidth',0.25); hold on;
        r_ci = abs(r_ci(1));
        if acos(dot(a,b)/(norm(a)*norm(b))) < 0                         % If the contacting arc is concave, set r_ci as negative
            r_ci = -r_ci;
        end
        % Use contact-rectange half-width to find a hertzian stiffness constant
        R_i = (10^-3)*r_ci*X.r_p/(r_ci+X.r_p); T_max = 32000;                      % Equivalent radius, max torque 
        F_max = T_max/(X.E*(10^-3)*(X.N+1));                               % Maximum expected pin force (for line fitting)
        b = sqrt(4*F_max*R_i/(pi*(10^-3)*X.w*X.E_red));                          % Half-width of contact rectangle
        r_ci_m = r_ci*(10^-3); r_p_m = X.r_p*(10^-3);
        ki = (10^3)*(r_ci_m-sqrt(r_ci_m^2-b^2)+r_p_m-sqrt(r_p_m^2-b^2))/F_max;        % Linearise the relationship between force and deformation
    end

    % ----------------------------------------------------------------------------------------------------------------------------------

    function ext_force_distribution(X)                                  % Calculate the external pin force distribution and total backlash
        % Force distribution will change as more pins engage. Consider each case.
        X.Fc = []; syms gamma_gl;                                          % gamma_g represts the total backlash inc. deformation
        for i = 1:(N+1)                                                 % Suppose i pins are in contact with the OHC
            displacement = X.gamma_c(1)-X.gamma_c([1:i])+gamma_gl;       % The vector of contact displacement distances   
            Fc = hertz_stiffness([1:i]).*displacement.*X.hat_Fc(:,1:i); % Calculate the stiffness constant for each contact point
            m_arm = (X.p_e-X.p_c(:,1:i));                               % The moment arm vector from p_e to each contact point
            torque = [0; 0; 0];
            for j = 1:i
                torque = torque+cross([m_arm(:,j); 0],[Fc(:,j); 0])     % Sum the moments of contact forces for i pins
            end
            eq1 = X.load/2 == torque(3);                                % CHECK ME: OHC experiences half of input force
            gamma_gl = solve(eq1,gamma_gl);                               % Solve for the total backlash given i contacting pins
            if gamma_gl < (X.gamma_c(i+1)-X.gamma_c(1))                  % Work out if the next pin would be in contact given gamma_g
                break                                                   % If the next pin would not be in contact, the supposed number
            end                                                         % ...of contacting pins == i, and gamma_g is physically correct
        end
        X.gamma_g = gamma_gl
        X.Fc = subs(X.Fc,gamma_gl);
    end

    % ----------------------------------------------------------------------------------------------------------------------------------

    function int_force_distribution(X)                                  % Calculate the load on the internal roller bearings
                                                                        % Assume that the force component at 90* to the moment arm
                                                                        % ...is proportional to the length of the moment arm
    end

    % ----------------------------------------------------------------------------------------------------------------------------------

    function out = trans(~,vec,shift,angle)                               % Linear tranformation of coordinate frames
        out = zeros(size(vec));
        for i = 1:length(vec(1,:))
            out(:,i) = [cos(angle),-sin(angle); sin(angle),cos(angle)]*vec(:,i)+shift;
        end
    end

    % ----------------------------------------------------------------------------------------------------------------------------------
    
    function x = base_eqn(X,dr_e,dr_p,t)                                   % The base function of an OHC
        Rr_e = X.r_e+dr_e; Rr_p = X.r_p+dr_p;                   
        ups = sqrt(X.E^2*(X.N+1)^2-2*X.E*Rr_e*(X.N+1)*cos(X.N*t)+Rr_e^2);  % upsilon is an offset factor from a hypocycloidal profile
        x = [Rr_e*(1-Rr_p./ups).*cos(t)-X.E*(1-Rr_p*(X.N+1)./ups).*cos((X.N+1)*t);
             Rr_e*(1-Rr_p./ups).*sin(t)-X.E*(1-Rr_p*(X.N+1)./ups).*sin((X.N+1)*t);];
    end

    function v_c = contact_speed(X,p_c)
        d = p_c-X.p_e; d_norm = norm(d);
        v_c_hat = [0,-1;1,0]*(d/d_norm);
        v_c = X.w*(X.p_e + v_c_hat*d_norm/X.N);
    end
end
end
% Find the projected backlash for given contact angle estimate   
function angle = proj_angle(X,i,c)
    % Get the contact point for the given contact angle estimate                                      
    p_c = frame_trans(X,base_eqn(X,X.dr_e,X.dr_p,c),X.ecc_pos,-X.t/X.N);                                             
    % Find the intersection point between the pin and the projected arc
    r1 = norm(p_c - X.ecc_pos); c1 = X.ecc_pos; 
    r2 = X.r_p; c2 = X.pin_pos(:,i);
    [x,y] = circcirc(c1(1),c1(2),r1,c2(1),c2(2),r2);
    if norm(p_c-[x(1);y(1)]) <= norm(p_c-[x(2);y(2)])
        p = [x(1);y(1)];
    else
        p = [x(2);y(2)];
    end                                                                
    X.cont_pin_pt(:,i) = p;
    % plot(p(1),p(2),'g*','MarkerSize',6); % TURN ON FOR POINT ESTIMATE VISUALS
    % Calculate the angle between the vectors that subtend the arc of interest
    v1 = p_c-X.ecc_pos; v2 = p-X.ecc_pos;
    angle = atan2(v1(1)*v2(2)-v1(2)*v2(1),v1(1)*v2(1)+v1(2)*v2(2));
    angle(isnan(angle)) = Inf;
    X.bkl_arc(i) = angle*r1;              
end
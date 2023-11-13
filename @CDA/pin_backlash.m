% For the current instant t (angle) calcaulate the pin backlash angles
function pin_backlash(X)
    % Initial guesses for closest contact points are at the peaks of OHC                                             
    X.cont_angle = pi/X.N + (0:1:X.N).*(2*pi/X.N);
    % Find the baclash angles for the inital guesses
    for i = 1:(X.N+1)
        X.bkl_angle(i) = X.proj_angle(i,X.cont_angle(i));
    end
    start_res = 0.01;    % Angular sampling resolution
    % For each pin recursively calculate the backlash angle (bilateral optim.)
    for i = 1:length(X.cont_angle)                               
        % Skip calculations for pins that will never contact
        if X.bkl_angle(i) < 0 || (X.bkl_angle(i) > 2*pi/(X.N+1))
            X.bkl_angle(i) = Inf;
            X.bkl_arc(i) = Inf;
        else
            X.contact_finder(i,start_res)
        end
    end
    angle = X.bkl_angle;
    arc_length = X.bkl_arc;
    % Filter negative values
    %X.bkl_angle(X.bkl_angle < 0) = Inf;
    % plot(X.cont_pin_pt(1,:),X.cont_pin_pt(2,:),'r.','MarkerSize',6)
end
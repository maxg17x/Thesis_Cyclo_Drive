% Find the closes contact point given the current estimates
function contact_finder(X,i,res)
    cont_pt_gear_frame = base_eqn(X,X.dr_e,X.dr_p,X.cont_angle(i));                              
    X.cont_pt(:,i) = frame_trans(X,cont_pt_gear_frame,X.ecc_pos,-X.t/X.N);
    set(X.h_cont_pt,'xdata',X.cont_pt(1,:),'ydata',X.cont_pt(2,:));
    % plot([X.cont_pin_pt(1,i) X.cont_pt(1,i)], [X.cont_pin_pt(2,i), X.cont_pt(2,i)],'b','Linewidth',2);
    res_thresh = 10^(-7);                                              % Recursively improve estimates till res_thresh is reached
    if res > res_thresh
        temp = X.cont_angle(i);                                    % Check if the left or right points offer improvements
        left_in = temp+res; right_in = temp-res; curr_out = X.bkl_angle(i);
        left_out = X.proj_angle(i,left_in); right_out = X.proj_angle(i,right_in);
        if left_out < curr_out
            X.cont_angle(i) = left_in;
            X.bkl_angle(i) = left_out;
            X.contact_finder(i,res);
        elseif right_out < curr_out
            X.cont_angle(i) = right_in;
            X.bkl_angle(i) = right_out;
            X.contact_finder(i,res);
        else
            X.contact_finder(i,0.5*res);                            % If the left/right estimates are worse lower search resolution
        end
    end
end
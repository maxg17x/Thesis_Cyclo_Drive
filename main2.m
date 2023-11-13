clear; clc; clf;
X = OHC3(36, 4.5, 1.5, 15, 0, 0.04, 0, 0, 0, 0, 0);
res = 500;

for i = 0:res
    pause(0.002);
    X.theta_shift(3*pi*i/res);
end
X.plot_backlash();s
% res = 384*4;
% rot = 1/16;
% % Arrays for recording
% bkl_angles = zeros(X.N+1, rot*res);
% for i = 1:rot*res
%     theta_shift(X, i*2*pi/res);
%     X.pin_backlash();
%     bkl_angles(1:(X.N+1), i) = X.angle_cont;
% end
% f2 = figure;
% bkl_angles(isinf(bkl_angles)) = NaN
% % for i = 1:(X.N+1)
% %     plot(1:rot*res,rad2deg(bkl_angles(i,:))); hold on;
% % end
% plot(1:rot*res,rad2deg(bkl_angles(i,:))); hold on;
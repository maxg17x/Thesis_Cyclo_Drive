
function plot_clearance();
    geo_param1 = [36,4.5,15,1.5, 0, 0.01];
    geo_param2 = [36,4.5,15,1.5, 0, 0.02];
    geo_param3 = [36,4.5,15,1.5, 0, 0.05];
    geo_param4 = [36,4.5,15,1.5, 0, 0.1];
    geo_param5 = [36,4.5,15,1.5, 0, 0.2];
    mat_param = [0.28,0.28,190*10^9,190*10^9,2000,2000,10,0.002];
    X1 = CDA(geo_param1, mat_param);
    X2 = CDA(geo_param2, mat_param);
    X3 = CDA(geo_param3, mat_param);
    X4 = CDA(geo_param4, mat_param);
    X5 = CDA(geo_param5, mat_param);
    figure(5);
    hp{1} = plot(1:(X1.N+1),X1.bkl_arc); hold on;
    hp{2} = plot(1:(X2.N+1),X2.bkl_arc); hold on;
    hp{3} = plot(1:(X3.N+1),X3.bkl_arc); hold on;
    hp{4} = plot(1:(X4.N+1),X4.bkl_arc); hold on;
    hp{5} = plot(1:(X5.N+1),X5.bkl_arc); hold on;
    plot(1:(X1.N+1),X1.bkl_arc,'k.'); hold on;
    plot(1:(X2.N+1),X2.bkl_arc,'k.'); hold on;
    plot(1:(X3.N+1),X3.bkl_arc,'k.'); hold on;
    plot(1:(X4.N+1),X4.bkl_arc,'k.'); hold on;
    plot(1:(X5.N+1),X5.bkl_arc,'k.'); hold on;
    xlim([1 7]); ylim([0 1]);
    xlabel('Pin Number'); ylabel('Clearance to Pin (mm)');
    hleg = legend([hp{1}(1);hp{2}(1);hp{3}(1);hp{4}(1);hp{5}(1)],'0.01 mm','0.02 mm', '0.05 mm', '0.1   mm', '0.2   mm');
    set(get(hleg,'Title'),'String','Offset Mod. (mm)')
end
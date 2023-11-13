function plot_eff_vs_ratio()
    eff_values = zeros(5, 21-8);
    offsets = [0.008, 0.01, 0.012, 0.014, 0.016];
    mat_params = [0.28,0.28,190*10^9,190*10^9,2000,2000,10,0.002];
    for j = 1:5
        for i = 8:21
            geo_params = [36,4.5,i,21/i, 0, offsets(j)];
            X = CDA(geo_params, mat_params);
            eff_values(j,i-7) = X.efficiency(50);
        end
    end
    figure(3);
    for k = 1:5
        hp{k} = plot(8:21,eff_values(k,:)); hold on;
    end
    xlim([8 21]); ylim([0.85 0.95]);
    xlabel('Reduction Ratio'); ylabel('Efficiency');
    hleg = legend([hp{1}(1);hp{2}(1);hp{3}(1);hp{4}(1);hp{5}(1)],'0.008 mm','0.010 mm', '0.012 mm', '0.014 mm', '0.016 mm');
    set(get(hleg,'Title'),'String','Offset Mod. (mm)')
end
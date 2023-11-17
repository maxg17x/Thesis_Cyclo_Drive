% Read in data files
clear; clc;
in_data = readtable('Hypersen_t1.xlsx');
out_data = readtable('datum_t1.xlsx');

% Average multiple torque values to the nearest second
in_info = [convertToSeconds(cellfun(@(s) s(12:22), in_data.Var2, 'UniformOutput', false)), abs(in_data.Var1)];
[uniqueTimes1, ~, ic] = unique(in_info(:, 1));
in_info_unique = [uniqueTimes1, accumarray(ic, in_info(:, 2), [], @mean)];

out_info = [convertToSeconds(out_data.Var7), out_data.Var2, out_data.Var3];
[uniqueTimes2, ~, ic] = unique(out_info(:, 1));
out_info_unique = [uniqueTimes2, accumarray(ic, out_info(:, 2), [], @mean), accumarray(ic, out_info(:, 3), [], @mean)];

% Trim early and late values from in_info_unique
early_offset = out_info_unique(1,1) - in_info_unique(1,1) + 1;
late_offset = (in_info_unique(end,1) - out_info_unique(end,1));
in_info_unique = in_info_unique(early_offset:(end-late_offset),:);

input_speed = 15*mean(out_data.Var5(166:end))*30/pi;

% Raw Data plots
figure(1);
data_span = 0:1:(length(in_info_unique(:,1))-1);
figure(1);
subplot(1,2,1);
plot(data_span,in_info_unique(:,2));
title('Raw Input Torque Readings');
xlabel('Time (s)');
ylabel('Input Torque (Nm)');
ylim([0 1.4])
subplot(2,2,2);
plot(data_span,out_info_unique(:,3));
title('Raw Output Torque Readings');
xlabel('Time (s)');
ylabel('Output Torque (Nm)');
ylim([0 18]);
sgtitle('Raw Torque Readings (Mapped and Adjusted)')

% Filter Data
in_torque_filtered = sgolayfilt(in_info_unique(:,2), 3, 75);
in_torque_filtered = in_torque_filtered(50:(end-50));
out_torque_filtered = sgolayfilt(out_info_unique(:,3), 3, 51);
out_torque_filtered = out_torque_filtered(50:(end-50));
out_target_torque = out_info_unique(:,2);
out_target_torque_trimmed = out_target_torque(50:(end-50));
data_span_filtered = 0:1:(length(in_torque_filtered)-1);


% Filtered Data Plots
figure(2);
subplot(1,2,1);
plot(data_span_filtered,in_torque_filtered);
title('Input Torque Readings');
xlabel('Time (s)');
ylabel('Input Torque (Nm)');
ylim([0 1.4])
subplot(2,2,2);
plot(data_span_filtered,out_torque_filtered); hold on;
plot(data_span_filtered, out_target_torque_trimmed);
title('Output Torque Readings');
xlabel('Time (s)');
ylabel('Output Torque (Nm)');
ylim([0 18]);
sgtitle('Filtered Torque Readings (Savitzky-Golay)')

% Efficency Plot
figure(3);
eff_real = out_torque_filtered./(15*in_torque_filtered);
eff_target = out_target_torque_trimmed./(15*in_torque_filtered);
eff_data = [out_target_torque_trimmed, 100*sgolayfilt(eff_real, 3, 121)];

plot(eff_data(:,1),eff_data(:,2),'k','LineWidth',1.25); hold on;
title('Filterd Efficiency Curve (99.8014 RPM avg.)');
xlabel('Load Torque (Nm)'); ylabel('Efficiency (%)');
ylim([0 100]); xlim([0 15]); grid on;
xticks(0:1:15); yticks(0:5:100);


function secondsSince12PM = convertToSeconds(timeStrings)
    secondsSince12PM = zeros(size(timeStrings)); % Preallocate for efficiency

    for i = 1:length(timeStrings)
        % Split the time string into its components
        timeParts = strsplit(timeStrings{i}, ':');

        % Convert each part to a number
        hour = str2double(timeParts{1});
        minute = str2double(timeParts{2});
        second = str2double(timeParts{3});
        millisecond = str2double(timeParts{4});

        % Convert hours, minutes, and seconds to total seconds
        % Add millisecond rounded to the nearest second
        totalSeconds = hour * 3600 + minute * 60 + second + round(millisecond / 1000);

        % Adjust for 12 PM
        secondsSince12PM(i) = totalSeconds - 12 * 3600;
    end
end

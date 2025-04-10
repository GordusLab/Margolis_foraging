function output = Fig2b(OUTPUT, M)

% Extract experimental reorientation data and time
data = OUTPUT.Model_Data;  % Reorientation events (rows: time, cols: individuals)
time = OUTPUT.Time;  % Time vector (frames)

% Compute moving average
[mean_rate, std_rate, initial_rates,rates] = moving_average(data, time);
time = time(end - length(mean_rate) + 1:end);

% histogram
num_trials = size(data,2); % Number of trials
measured_values = rates; 
bin_edges = 0:0.5:3; % Bin edges
bin_centers = (bin_edges(1:end-1) + bin_edges(2:end)) / 2;
bin_num = length(bin_edges)-1;
% Compute Histogram Counts (for each time point)
hist_counts = zeros(length(time),bin_num);
for i = 1:length(time)
    hist_counts(i, :) = histcounts(measured_values(i, :), bin_edges);
end

% Normalize Counts to [0, 1] for Colormap
norm_counts = hist_counts / max(hist_counts(:)); % Normalize

% Plot
imagesc(time, bin_centers, norm_counts'); % Transpose for correct orientation
axis xy; % Flip the y-axis for proper alignment
colormap([linspace(1, 1, 256)', linspace(1, 0, 256)', linspace(1, 0, 256)']); % White to red
colorbar;
xlabel('Time');
ylabel('Measured Value');
xlim([0, max(time)]);
title('2D Histogram of Measurements');

% Plotting mean_rate
% Extract relevant data
hold on;
plot(time, mean_rate, 'b', 'LineWidth', 1.5); hold on;
xlabel('Time (minutes)');
ylabel('Average Reorientation Rate');
title('Moving Average of Reorientation Rate Over Time');
grid on;


% Extract relevant data
hold on;
plot(time, mean_rate, 'k', 'LineWidth', 1.5); hold on;
xlabel('Time (minutes)');
ylabel('Average Reorientation Rate');
title('Moving Average of Reorientation Rate Over Time');
grid on;

% Generate filename dynamically: e.g., "M1000_Trial1_Model_Output.mat"
filename = sprintf('M%d_Fig2b_lower.fig', M);
foldername = sprintf('M%d', M);
filename = fullfile(strcat('Figures/',foldername), filename);
saveas(gcf, filename);
close(gcf);

% Extract relevant data
figure;
plot(time, mean_rate, 'r', 'LineWidth', 1.5); hold on;
plot(time, mean_rate+std_rate, 'r--', 'LineWidth', 1.5); hold on;
plot(time, mean_rate-std_rate, 'r--', 'LineWidth', 1.5); hold on;

xlabel('Time (minutes)');
ylabel('Average Reorientation Rate');
title('Moving Average of Reorientation Rate Over Time');
grid on;

data = OUTPUT.Exp_Data;  % Reorientation events (rows: time, cols: individuals)
time = OUTPUT.Time;  % Time vector (frames)

% EXPERIMENTAL DATA

% Compute moving average
[mean_rate, std_rate, initial_rates,rates] = moving_average(data, time);
time = time(end - length(mean_rate) + 1:end);

hold on;
plot(time, mean_rate, 'b', 'LineWidth', 1.5); hold on;
plot(time, mean_rate+std_rate, 'b--', 'LineWidth', 1.5); hold on;
plot(time, mean_rate-std_rate, 'b--', 'LineWidth', 1.5); hold on;

xlabel('Time (minutes)');
ylabel('Average Reorientation Rate');
title('Moving Average of Reorientation Rate Over Time');
grid on;

filename = sprintf('M%d_Fig2b_upper.fig', M);
foldername = sprintf('M%d', M);
filename = fullfile(strcat('Figures/',foldername), filename);
saveas(gcf, filename);
close(gcf);
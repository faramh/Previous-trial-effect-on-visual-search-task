% PSTH Approach for vlPFC
%% Load Results of separation
clear all; close all; 
% Load results from each part
load('results_part1.mat');
load('results_part2.mat');
load('results_part3.mat');

% Combine results from all parts
TP_prev_TP_Final = [TP_prev_TP; TP_prev_TP_Part2; TP_prev_TP_Part3];
TP_prev_TP_Pre_Final = [TP_prev_TP_Pre; TP_prev_TP_Pre_Part2; TP_prev_TP_Pre_Part3];
TP_prev_TA_Final = [TP_prev_TA; TP_prev_TA_Part2; TP_prev_TA_Part3];
TP_prev_TA_Pre_Final = [TP_prev_TA_Pre; TP_prev_TA_Pre_Part2; TP_prev_TA_Pre_Part3];
TA_prev_TP_Final = [TA_prev_TP; TA_prev_TP_Part2; TA_prev_TP_Part3];
TA_prev_TP_Pre_Final = [TA_prev_TP_Pre; TA_prev_TP_Pre_Part2; TA_prev_TP_Pre_Part3];
TA_prev_TA_Final = [TA_prev_TA; TA_prev_TA_Part2; TA_prev_TA_Part3];
TA_prev_TA_Pre_Final = [TA_prev_TA_Pre; TA_prev_TA_Pre_Part2; TA_prev_TA_Pre_Part3];


%% PSTH vlPFC
 

column_names = cell(1, 1600);

for i = 1:1600
    column_names{i} = ['bin', num2str(i)];
end



% TP-TP
TP_prev_TP_Final_bins = TP_prev_TP_Final{:, column_names(1:1600)};
TP_prev_TP_Final_bins_mean=nanmean(TP_prev_TP_Final_bins,1);
% TP-TA
TP_prev_TA_Final_bins = TP_prev_TA_Final{:, column_names(1:1600)};
TP_prev_TA_Final_bins_mean=nanmean(TP_prev_TA_Final_bins,1);
%TA-TP
TA_prev_TP_Final_bins = TA_prev_TP_Final{:, column_names(1:1600)};
TA_prev_TP_Final_bins_mean=nanmean(TA_prev_TP_Final_bins,1);
%TA-TA
TA_prev_TA_Final_bins = TA_prev_TA_Final{:, column_names(1:1600)};
TA_prev_TA_Final_bins_mean=nanmean(TA_prev_TA_Final_bins,1);


% Define time axis (assuming your data spans 1 second)
time_axis = linspace(-0.6, 1, 1600);

% Plot firing rates for each condition
figure;
hold on;
colors = get_distinguishable_colors(8); % You can also use other colormap functions like parula, jet, etc.
 
PlotPSTH(time_axis,TP_prev_TP_Final_bins_mean,colors(1,:));
PlotPSTH(time_axis,TP_prev_TA_Final_bins_mean,colors(2,:));

% Labeling
xlabel('Time (s)');
ylabel('Firing Rate (spikes/s)');
title('PSTH Firing Rates Across Conditions');
legend('TP Prev TP','TP prev TA');
grid on;
hold off;


figure;
hold on;
colors = get_distinguishable_colors(8); % You can also use other colormap functions like parula, jet, etc.
 
PlotPSTH(time_axis,TA_prev_TP_Final_bins_mean,colors(1,:));
PlotPSTH(time_axis,TA_prev_TA_Final_bins_mean,colors(2,:));

% Labeling
xlabel('Time (s)');
ylabel('Firing Rate (spikes/s)');
title('PSTH Firing Rates Across Conditions');
legend('TA Prev TP','TA prev TA');

grid on;
hold off;




%% Separate monkey


% Get data for Monkey 1
Monkey1_TP_prev_TP_Final = TP_prev_TP_Final(TP_prev_TP_Final.Subject == 1, :);
Monkey1_TP_prev_TA_Final = TP_prev_TA_Final(TP_prev_TA_Final.Subject == 1, :);
Monkey1_TA_prev_TP_Final = TA_prev_TP_Final(TA_prev_TP_Final.Subject == 1, :);
Monkey1_TA_prev_TA_Final = TA_prev_TA_Final(TA_prev_TA_Final.Subject == 1, :);

% Get data for Monkey 2
Monkey2_TP_prev_TP_Final = TP_prev_TP_Final(TP_prev_TP_Final.Subject == 2, :);
Monkey2_TP_prev_TA_Final = TP_prev_TA_Final(TP_prev_TA_Final.Subject == 2, :);
Monkey2_TA_prev_TP_Final = TA_prev_TP_Final(TA_prev_TP_Final.Subject == 2, :);
Monkey2_TA_prev_TA_Final = TA_prev_TA_Final(TA_prev_TA_Final.Subject == 2, :);

% Define the columns (bins)
column_names = cell(1, 1600);
for i = 1:1600
    column_names{i} = ['bin', num2str(i)];
end

% Compute mean firing rates for Monkey 1
Monkey1_TP_prev_TP_Final_bins = Monkey1_TP_prev_TP_Final{:, column_names(1:1600)};
Monkey1_TP_prev_TP_Final_bins_mean = nanmean(Monkey1_TP_prev_TP_Final_bins, 1);

Monkey1_TP_prev_TA_Final_bins = Monkey1_TP_prev_TA_Final{:, column_names(1:1600)};
Monkey1_TP_prev_TA_Final_bins_mean = nanmean(Monkey1_TP_prev_TA_Final_bins, 1);

Monkey1_TA_prev_TP_Final_bins = Monkey1_TA_prev_TP_Final{:, column_names(1:1600)};
Monkey1_TA_prev_TP_Final_bins_mean = nanmean(Monkey1_TA_prev_TP_Final_bins, 1);

Monkey1_TA_prev_TA_Final_bins = Monkey1_TA_prev_TA_Final{:, column_names(1:1600)};
Monkey1_TA_prev_TA_Final_bins_mean = nanmean(Monkey1_TA_prev_TA_Final_bins, 1);

% Define the time axis
time_axis = linspace(-0.6, 1, 1600);

% Plot Monkey 1 firing rates
figure;
hold on;
colors = get_distinguishable_colors(8);
PlotPSTH(time_axis, Monkey1_TP_prev_TP_Final_bins_mean, colors(1,:));
PlotPSTH(time_axis, Monkey1_TP_prev_TA_Final_bins_mean, colors(2,:));

xlabel('Time (s)');
ylabel('Firing Rate (spikes/s)');
title('Monkey 1: PSTH Firing Rates Across Conditions');
legend('TP Prev TP', 'TP Prev TA');
grid on;
hold off;

figure;
hold on;
PlotPSTH(time_axis, Monkey1_TA_prev_TP_Final_bins_mean, colors(1,:));
PlotPSTH(time_axis, Monkey1_TA_prev_TA_Final_bins_mean, colors(2,:));

xlabel('Time (s)');
ylabel('Firing Rate (spikes/s)');
title('Monkey 1: PSTH Firing Rates Across Conditions');
legend('TA Prev TP', 'TA Prev TA');
grid on;
hold off;

% Repeat the same for Monkey 2
Monkey2_TP_prev_TP_Final_bins = Monkey2_TP_prev_TP_Final{:, column_names(1:1600)};
Monkey2_TP_prev_TP_Final_bins_mean = nanmean(Monkey2_TP_prev_TP_Final_bins, 1);

Monkey2_TP_prev_TA_Final_bins = Monkey2_TP_prev_TA_Final{:, column_names(1:1600)};
Monkey2_TP_prev_TA_Final_bins_mean = nanmean(Monkey2_TP_prev_TA_Final_bins, 1);

Monkey2_TA_prev_TP_Final_bins = Monkey2_TA_prev_TP_Final{:, column_names(1:1600)};
Monkey2_TA_prev_TP_Final_bins_mean = nanmean(Monkey2_TA_prev_TP_Final_bins, 1);

Monkey2_TA_prev_TA_Final_bins = Monkey2_TA_prev_TA_Final{:, column_names(1:1600)};
Monkey2_TA_prev_TA_Final_bins_mean = nanmean(Monkey2_TA_prev_TA_Final_bins, 1);

% Plot Monkey 2 firing rates
figure;
hold on;
PlotPSTH(time_axis, Monkey2_TP_prev_TP_Final_bins_mean, colors(1,:));
PlotPSTH(time_axis, Monkey2_TP_prev_TA_Final_bins_mean, colors(2,:));

xlabel('Time (s)');
ylabel('Firing Rate (spikes/s)');
title('Monkey 2: PSTH Firing Rates Across Conditions');
legend('TP Prev TP', 'TP Prev TA');
grid on;
hold off;

figure;
hold on;
PlotPSTH(time_axis, Monkey2_TA_prev_TP_Final_bins_mean, colors(1,:));
PlotPSTH(time_axis, Monkey2_TA_prev_TA_Final_bins_mean, colors(2,:));

xlabel('Time (s)');
ylabel('Firing Rate (spikes/s)');
title('Monkey 2: PSTH Firing Rates Across Conditions');
legend('TA Prev TP', 'TA Prev TA');
grid on;
hold off;





%% % Separate Eff and Ineff conditions 
% Eff_TPTA = find(TP_prev_TP_Final.SlopeTP<=20);
% Ineff_TPTA = find(TP_prev_TP_Final.SlopeTP>=35);

% Separate Eff and Ineff conditions for TA_prev_TP_Final
%TP - TP
Eff_TPTP = find(TP_prev_TP_Final.SlopeTP<=20);
Ineff_TPTP = find(TP_prev_TP_Final.SlopeTP>=35);


TP_prev_TP_Final_Eff=TP_prev_TP_Final(Eff_TPTP,:);
TP_prev_TP_Final_Ineff=TP_prev_TP_Final(Ineff_TPTP,:);

%TP - TA
Eff_TPTA = find(TP_prev_TA_Final.SlopeTP<=20);
Ineff_TPTA = find(TP_prev_TA_Final.SlopeTP>=35);

TP_prev_TA_Final_Eff=TP_prev_TA_Final(Eff_TPTA,:);
TP_prev_TA_Final_Ineff=TP_prev_TA_Final(Ineff_TPTA,:);

%TA - TP
Eff_TATP = find(TA_prev_TP_Final.SlopeTP<=20);
Ineff_TATP = find(TA_prev_TP_Final.SlopeTP>=35);




TA_prev_TP_Final_Eff=TA_prev_TP_Final(Eff_TATP,:);
TA_prev_TP_Final_Ineff=TA_prev_TP_Final(Ineff_TATP,:);

%TA - TA
Eff_TATA = find(TA_prev_TA_Final.SlopeTP<=20);
Ineff_TATA = find(TA_prev_TA_Final.SlopeTP>=35);

TA_prev_TA_Final_Eff=TA_prev_TA_Final(Eff_TATA,:);
TA_prev_TA_Final_Ineff=TA_prev_TA_Final(Ineff_TATA,:);


% Accessing to bins of PSTH
% TP-TP Eff
TP_prev_TP_Final_Eff_bins = TP_prev_TP_Final_Eff{:, column_names(1:1600)};
TP_prev_TP_Final_Eff_bins_mean = nanmean(TP_prev_TP_Final_Eff_bins, 1);

% TP-TP Ineff
TP_prev_TP_Final_Ineff_bins = TP_prev_TP_Final_Ineff{:, column_names(1:1600)};
TP_prev_TP_Final_Ineff_bins_mean = nanmean(TP_prev_TP_Final_Ineff_bins, 1);

% TP-TA Eff
TP_prev_TA_Final_Eff_bins = TP_prev_TA_Final_Eff{:, column_names(1:1600)};
TP_prev_TA_Final_Eff_bins_mean = nanmean(TP_prev_TA_Final_Eff_bins, 1);

% TP-TA Ineff
TP_prev_TA_Final_Ineff_bins = TP_prev_TA_Final_Ineff{:, column_names(1:1600)};
TP_prev_TA_Final_Ineff_bins_mean = nanmean(TP_prev_TA_Final_Ineff_bins, 1);

% TA-TP Eff
TA_prev_TP_Final_Eff_bins = TA_prev_TP_Final_Eff{:, column_names(1:1600)};
TA_prev_TP_Final_Eff_bins_mean = nanmean(TA_prev_TP_Final_Eff_bins, 1);

% TA-TP Ineff
TA_prev_TP_Final_Ineff_bins = TA_prev_TP_Final_Ineff{:, column_names(1:1600)};
TA_prev_TP_Final_Ineff_bins_mean = nanmean(TA_prev_TP_Final_Ineff_bins, 1);

% TA-TA Eff
TA_prev_TA_Final_Eff_bins = TA_prev_TA_Final_Eff{:, column_names(1:1600)};
TA_prev_TA_Final_Eff_bins_mean = nanmean(TA_prev_TA_Final_Eff_bins, 1);

% TA-TA Ineff
TA_prev_TA_Final_Ineff_bins = TA_prev_TA_Final_Ineff{:, column_names(1:1600)};
TA_prev_TA_Final_Ineff_bins_mean = nanmean(TA_prev_TA_Final_Ineff_bins, 1);





figure;
hold on;
colors = get_distinguishable_colors(8); % You can also use other colormap functions like parula, jet, etc.
 
PlotPSTH(time_axis,TP_prev_TP_Final_Eff_bins_mean,colors(1,:));
PlotPSTH(time_axis,TP_prev_TA_Final_Eff_bins_mean,colors(2,:));

% Labeling
xlabel('Time (s)');
ylabel('Firing Rate (spikes/s)');
title('PSTH Firing Rates Across Conditions for Effiecient');
legend('TP Prev TP','TP prev TA');
grid on;
hold off;


figure;
hold on;
colors = get_distinguishable_colors(8); % You can also use other colormap functions like parula, jet, etc.
 
PlotPSTH(time_axis,TP_prev_TP_Final_Ineff_bins_mean,colors(1,:));
PlotPSTH(time_axis,TP_prev_TA_Final_Ineff_bins_mean,colors(2,:));

% Labeling
xlabel('Time (s)');
ylabel('Firing Rate (spikes/s)');
title('PSTH Firing Rates Across Conditions for Inefficient');
legend('TP Prev TP','TP prev TA');

grid on;
hold off;

figure;
hold on;
colors = get_distinguishable_colors(8); % You can also use other colormap functions like parula, jet, etc.
 
PlotPSTH(time_axis,TA_prev_TP_Final_Eff_bins_mean,colors(1,:));
PlotPSTH(time_axis,TA_prev_TA_Final_Eff_bins_mean,colors(2,:));

% Labeling
xlabel('Time (s)');
ylabel('Firing Rate (spikes/s)');
title('PSTH Firing Rates Across Conditions Efficient');
legend('TA Prev TP','TA prev TA');
grid on;
hold off;


figure;
hold on;
colors = get_distinguishable_colors(8); % You can also use other colormap functions like parula, jet, etc.
 
PlotPSTH(time_axis,TA_prev_TP_Final_Ineff_bins_mean,colors(1,:));
PlotPSTH(time_axis,TA_prev_TA_Final_Ineff_bins_mean,colors(2,:));

% Labeling
xlabel('Time (s)');
ylabel('Firing Rate (spikes/s)');
title('PSTH Firing Rates Across Conditions for Inefficient');
legend('TA Prev TP','TA prev TA');

grid on;
hold off;

%% Monkey sep Eff/INeff

% Separate data by Subject (Monkey 1 and Monkey 2)
% Monkey 1
Monkey1_TP_prev_TP_Final = TP_prev_TP_Final(TP_prev_TP_Final.Subject == 1, :);
Monkey1_TP_prev_TA_Final = TP_prev_TA_Final(TP_prev_TA_Final.Subject == 1, :);
Monkey1_TA_prev_TP_Final = TA_prev_TP_Final(TA_prev_TP_Final.Subject == 1, :);
Monkey1_TA_prev_TA_Final = TA_prev_TA_Final(TA_prev_TA_Final.Subject == 1, :);

% Monkey 2
Monkey2_TP_prev_TP_Final = TP_prev_TP_Final(TP_prev_TP_Final.Subject == 2, :);
Monkey2_TP_prev_TA_Final = TP_prev_TA_Final(TP_prev_TA_Final.Subject == 2, :);
Monkey2_TA_prev_TP_Final = TA_prev_TP_Final(TA_prev_TP_Final.Subject == 2, :);
Monkey2_TA_prev_TA_Final = TA_prev_TA_Final(TA_prev_TA_Final.Subject == 2, :);

% For each monkey, separate Eff and Ineff conditions
% Monkey 1
Eff_TPTP_Monkey1 = find(Monkey1_TP_prev_TP_Final.SlopeTP <= 20);
Ineff_TPTP_Monkey1 = find(Monkey1_TP_prev_TP_Final.SlopeTP >= 35);
Monkey1_TP_prev_TP_Final_Eff = Monkey1_TP_prev_TP_Final(Eff_TPTP_Monkey1, :);
Monkey1_TP_prev_TP_Final_Ineff = Monkey1_TP_prev_TP_Final(Ineff_TPTP_Monkey1, :);

Eff_TPTA_Monkey1 = find(Monkey1_TP_prev_TA_Final.SlopeTP <= 20);
Ineff_TPTA_Monkey1 = find(Monkey1_TP_prev_TA_Final.SlopeTP >= 35);
Monkey1_TP_prev_TA_Final_Eff = Monkey1_TP_prev_TA_Final(Eff_TPTA_Monkey1, :);
Monkey1_TP_prev_TA_Final_Ineff = Monkey1_TP_prev_TA_Final(Ineff_TPTA_Monkey1, :);

Eff_TATP_Monkey1 = find(Monkey1_TA_prev_TP_Final.SlopeTP <= 20);
Ineff_TATP_Monkey1 = find(Monkey1_TA_prev_TP_Final.SlopeTP >= 35);
Monkey1_TA_prev_TP_Final_Eff = Monkey1_TA_prev_TP_Final(Eff_TATP_Monkey1, :);
Monkey1_TA_prev_TP_Final_Ineff = Monkey1_TA_prev_TP_Final(Ineff_TATP_Monkey1, :);

Eff_TATA_Monkey1 = find(Monkey1_TA_prev_TA_Final.SlopeTP <= 20);
Ineff_TATA_Monkey1 = find(Monkey1_TA_prev_TA_Final.SlopeTP >= 35);
Monkey1_TA_prev_TA_Final_Eff = Monkey1_TA_prev_TA_Final(Eff_TATA_Monkey1, :);
Monkey1_TA_prev_TA_Final_Ineff = Monkey1_TA_prev_TA_Final(Ineff_TATA_Monkey1, :);

% Monkey 2
Eff_TPTP_Monkey2 = find(Monkey2_TP_prev_TP_Final.SlopeTP <= 20);
Ineff_TPTP_Monkey2 = find(Monkey2_TP_prev_TP_Final.SlopeTP >= 35);
Monkey2_TP_prev_TP_Final_Eff = Monkey2_TP_prev_TP_Final(Eff_TPTP_Monkey2, :);
Monkey2_TP_prev_TP_Final_Ineff = Monkey2_TP_prev_TP_Final(Ineff_TPTP_Monkey2, :);

Eff_TPTA_Monkey2 = find(Monkey2_TP_prev_TA_Final.SlopeTP <= 20);
Ineff_TPTA_Monkey2 = find(Monkey2_TP_prev_TA_Final.SlopeTP >= 35);
Monkey2_TP_prev_TA_Final_Eff = Monkey2_TP_prev_TA_Final(Eff_TPTA_Monkey2, :);
Monkey2_TP_prev_TA_Final_Ineff = Monkey2_TP_prev_TA_Final(Ineff_TPTA_Monkey2, :);

Eff_TATP_Monkey2 = find(Monkey2_TA_prev_TP_Final.SlopeTP <= 20);
Ineff_TATP_Monkey2 = find(Monkey2_TA_prev_TP_Final.SlopeTP >= 35);
Monkey2_TA_prev_TP_Final_Eff = Monkey2_TA_prev_TP_Final(Eff_TATP_Monkey2, :);
Monkey2_TA_prev_TP_Final_Ineff = Monkey2_TA_prev_TP_Final(Ineff_TATP_Monkey2, :);

Eff_TATA_Monkey2 = find(Monkey2_TA_prev_TA_Final.SlopeTP <= 20);
Ineff_TATA_Monkey2 = find(Monkey2_TA_prev_TA_Final.SlopeTP >= 35);
Monkey2_TA_prev_TA_Final_Eff = Monkey2_TA_prev_TA_Final(Eff_TATA_Monkey2, :);
Monkey2_TA_prev_TA_Final_Ineff = Monkey2_TA_prev_TA_Final(Ineff_TATA_Monkey2, :);

%% Plot 
% Accessing bins of PSTH for Monkey 1
% TP-TP Eff Monkey 1
Monkey1_TP_prev_TP_Final_Eff_bins = Monkey1_TP_prev_TP_Final_Eff{:, column_names(1:1600)};
Monkey1_TP_prev_TP_Final_Eff_bins_mean = nanmean(Monkey1_TP_prev_TP_Final_Eff_bins, 1);

% TP-TP Ineff Monkey 1
Monkey1_TP_prev_TP_Final_Ineff_bins = Monkey1_TP_prev_TP_Final_Ineff{:, column_names(1:1600)};
Monkey1_TP_prev_TP_Final_Ineff_bins_mean = nanmean(Monkey1_TP_prev_TP_Final_Ineff_bins, 1);

% TP-TA Eff Monkey 1
Monkey1_TP_prev_TA_Final_Eff_bins = Monkey1_TP_prev_TA_Final_Eff{:, column_names(1:1600)};
Monkey1_TP_prev_TA_Final_Eff_bins_mean = nanmean(Monkey1_TP_prev_TA_Final_Eff_bins, 1);

% TP-TA Ineff Monkey 1
Monkey1_TP_prev_TA_Final_Ineff_bins = Monkey1_TP_prev_TA_Final_Ineff{:, column_names(1:1600)};
Monkey1_TP_prev_TA_Final_Ineff_bins_mean = nanmean(Monkey1_TP_prev_TA_Final_Ineff_bins, 1);

% TA-TP Eff Monkey 1
Monkey1_TA_prev_TP_Final_Eff_bins = Monkey1_TA_prev_TP_Final_Eff{:, column_names(1:1600)};
Monkey1_TA_prev_TP_Final_Eff_bins_mean = nanmean(Monkey1_TA_prev_TP_Final_Eff_bins, 1);

% TA-TP Ineff Monkey 1
Monkey1_TA_prev_TP_Final_Ineff_bins = Monkey1_TA_prev_TP_Final_Ineff{:, column_names(1:1600)};
Monkey1_TA_prev_TP_Final_Ineff_bins_mean = nanmean(Monkey1_TA_prev_TP_Final_Ineff_bins, 1);

% TA-TA Eff Monkey 1
Monkey1_TA_prev_TA_Final_Eff_bins = Monkey1_TA_prev_TA_Final_Eff{:, column_names(1:1600)};
Monkey1_TA_prev_TA_Final_Eff_bins_mean = nanmean(Monkey1_TA_prev_TA_Final_Eff_bins, 1);

% TA-TA Ineff Monkey 1
Monkey1_TA_prev_TA_Final_Ineff_bins = Monkey1_TA_prev_TA_Final_Ineff{:, column_names(1:1600)};
Monkey1_TA_prev_TA_Final_Ineff_bins_mean = nanmean(Monkey1_TA_prev_TA_Final_Ineff_bins, 1);

% Now for Monkey 2
% TP-TP Eff Monkey 2
Monkey2_TP_prev_TP_Final_Eff_bins = Monkey2_TP_prev_TP_Final_Eff{:, column_names(1:1600)};
Monkey2_TP_prev_TP_Final_Eff_bins_mean = nanmean(Monkey2_TP_prev_TP_Final_Eff_bins, 1);

% TP-TP Ineff Monkey 2
Monkey2_TP_prev_TP_Final_Ineff_bins = Monkey2_TP_prev_TP_Final_Ineff{:, column_names(1:1600)};
Monkey2_TP_prev_TP_Final_Ineff_bins_mean = nanmean(Monkey2_TP_prev_TP_Final_Ineff_bins, 1);

% TP-TA Eff Monkey 2
Monkey2_TP_prev_TA_Final_Eff_bins = Monkey2_TP_prev_TA_Final_Eff{:, column_names(1:1600)};
Monkey2_TP_prev_TA_Final_Eff_bins_mean = nanmean(Monkey2_TP_prev_TA_Final_Eff_bins, 1);

% TP-TA Ineff Monkey 2
Monkey2_TP_prev_TA_Final_Ineff_bins = Monkey2_TP_prev_TA_Final_Ineff{:, column_names(1:1600)};
Monkey2_TP_prev_TA_Final_Ineff_bins_mean = nanmean(Monkey2_TP_prev_TA_Final_Ineff_bins, 1);

% TA-TP Eff Monkey 2
Monkey2_TA_prev_TP_Final_Eff_bins = Monkey2_TA_prev_TP_Final_Eff{:, column_names(1:1600)};
Monkey2_TA_prev_TP_Final_Eff_bins_mean = nanmean(Monkey2_TA_prev_TP_Final_Eff_bins, 1);

% TA-TP Ineff Monkey 2
Monkey2_TA_prev_TP_Final_Ineff_bins = Monkey2_TA_prev_TP_Final_Ineff{:, column_names(1:1600)};
Monkey2_TA_prev_TP_Final_Ineff_bins_mean = nanmean(Monkey2_TA_prev_TP_Final_Ineff_bins, 1);

% TA-TA Eff Monkey 2
Monkey2_TA_prev_TA_Final_Eff_bins = Monkey2_TA_prev_TA_Final_Eff{:, column_names(1:1600)};
Monkey2_TA_prev_TA_Final_Eff_bins_mean = nanmean(Monkey2_TA_prev_TA_Final_Eff_bins, 1);

% TA-TA Ineff Monkey 2
Monkey2_TA_prev_TA_Final_Ineff_bins = Monkey2_TA_prev_TA_Final_Ineff{:, column_names(1:1600)};
Monkey2_TA_prev_TA_Final_Ineff_bins_mean = nanmean(Monkey2_TA_prev_TA_Final_Ineff_bins, 1);

% Plot for Monkey 1
figure;
hold on;
colors = get_distinguishable_colors(8);
PlotPSTH(time_axis, Monkey1_TP_prev_TP_Final_Eff_bins_mean, colors(1,:));
PlotPSTH(time_axis, Monkey1_TP_prev_TA_Final_Eff_bins_mean, colors(2,:));
xlabel('Time (s)');
ylabel('Firing Rate (spikes/s)');
title('Monkey 1: PSTH Firing Rates Across Conditions for Efficient');
legend('TP Prev TP','TP Prev TA');
grid on;
hold off;

figure;
hold on;
PlotPSTH(time_axis, Monkey1_TP_prev_TP_Final_Ineff_bins_mean, colors(1,:));
PlotPSTH(time_axis, Monkey1_TP_prev_TA_Final_Ineff_bins_mean, colors(2,:));
xlabel('Time (s)');
ylabel('Firing Rate (spikes/s)');
title('Monkey 1: PSTH Firing Rates Across Conditions for Inefficient');
legend('TP Prev TP','TP Prev TA');
grid on;
hold off;

figure;
hold on;
PlotPSTH(time_axis, Monkey1_TA_prev_TP_Final_Eff_bins_mean, colors(1,:));
PlotPSTH(time_axis, Monkey1_TA_prev_TA_Final_Eff_bins_mean, colors(2,:));
xlabel('Time (s)');
ylabel('Firing Rate (spikes/s)');
title('Monkey 1: PSTH Firing Rates Across Conditions for Efficient');
legend('TA Prev TP','TA Prev TA');
grid on;
hold off;

figure;
hold on;
PlotPSTH(time_axis, Monkey1_TA_prev_TP_Final_Ineff_bins_mean, colors(1,:));
PlotPSTH(time_axis, Monkey1_TA_prev_TA_Final_Ineff_bins_mean, colors(2,:));
xlabel('Time (s)');
ylabel('Firing Rate (spikes/s)');
title('Monkey 1: PSTH Firing Rates Across Conditions for Inefficient');
legend('TA Prev TP','TA Prev TA');
grid on;
hold off;

% Plot for Monkey 2
figure;
hold on;
PlotPSTH(time_axis, Monkey2_TP_prev_TP_Final_Eff_bins_mean, colors(1,:));
PlotPSTH(time_axis, Monkey2_TP_prev_TA_Final_Eff_bins_mean, colors(2,:));
xlabel('Time (s)');
ylabel('Firing Rate (spikes/s)');
title('Monkey 2: PSTH Firing Rates Across Conditions for Efficient');
legend('TP Prev TP','TP Prev TA');
grid on;
hold off;

figure;
hold on;
PlotPSTH(time_axis, Monkey2_TP_prev_TP_Final_Ineff_bins_mean, colors(1,:));
PlotPSTH(time_axis, Monkey2_TP_prev_TA_Final_Ineff_bins_mean, colors(2,:));
xlabel('Time (s)');
ylabel('Firing Rate (spikes/s)');
title('Monkey 2: PSTH Firing Rates Across Conditions for Inefficient');
legend('TP Prev TP','TP Prev TA');
grid on;
hold off;

figure;
hold on;
PlotPSTH(time_axis, Monkey2_TA_prev_TP_Final_Eff_bins_mean, colors(1,:));
PlotPSTH(time_axis, Monkey2_TA_prev_TA_Final_Eff_bins_mean, colors(2,:));
xlabel('Time (s)');
ylabel('Firing Rate (spikes/s)');
title('Monkey 2: PSTH Firing Rates Across Conditions for Efficient');
legend('TA Prev TP','TA Prev TA');
grid on;
hold off;

figure;
hold on;
PlotPSTH(time_axis, Monkey2_TA_prev_TP_Final_Ineff_bins_mean, colors(1,:));
PlotPSTH(time_axis, Monkey2_TA_prev_TA_Final_Ineff_bins_mean, colors(2,:));
xlabel('Time (s)');
ylabel('Firing Rate (spikes/s)');
title('Monkey 2: PSTH Firing Rates Across Conditions for Inefficient');
legend('TA Prev TP','TA Prev TA');
grid on;
hold off;






%% Functions





% Plot PSTH
% Gives mean bins of the data 
function PlotPSTH(time_axis,meanbins,color)

sigma1 = 4; % Adjust sigma value as needed
 % Adjust window size as needed
smoothed_data = imgaussfilt(meanbins, sigma1);

% Plot PSTH as a continuous line

plot(time_axis, smoothed_data, 'Color', color, 'LineWidth', 2);

end


function colors = get_distinguishable_colors(n)
    % Generate a set of distinguishable colors
    golden_ratio_conjugate = (1 + sqrt(5)) / 2;
    hue = mod((0:n-1)' / golden_ratio_conjugate, 1);
    saturation = 0.6;
    value = 0.9;
    colors = hsv2rgb([hue, saturation * ones(n, 1), value * ones(n, 1)]);
end




















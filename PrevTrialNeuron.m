% In this part of analyze we want to consider Behavior
close all; clear all;
addpath("/Users/fara/Desktop/Breath Deeper/Current Projects/Previous Trial Effect/Neuron Analysis")
% SNr Data - Comment/Uncomment
load("SNr_new.mat");
data=table;
% vlPFC Data - Uncomment/Comment
% load("vlPFC.mat"); 
% data=tableB;

%% Condition separation
% C= Current P=Previous
% Conditions: C:TP/P:TP -- C:TP/P:TA -- C:TA/P:TP -- C:TA/P:TA


% Initialize arrays to hold the concatenated results
TP_prev_TP = [];
TP_prev_TA = [];
TA_prev_TP = [];
TA_prev_TA = [];

TP_prev_TP_Pre = [];
TP_prev_TA_Pre = [];
TA_prev_TP_Pre = [];
TA_prev_TA_Pre = [];

counter = 0;

% Loop through the data to detect session changes based on TaskOrder
prevTaskOrder = data.TaskOrder(1); % Initialize with the first TaskOrder
session_start_idx = 1; % Start index for the first session
Valids=0;
for i = 2:height(data)
    currentTaskOrder = data.TaskOrder(i);
    disp(i);
    
    % Detect a change in TaskOrder, indicating a new session
    if currentTaskOrder ~= prevTaskOrder
        % Extract session data from session_start_idx to i-1
        session_data = data(session_start_idx:i-1, :);
        TrialsIndices = session_data.TrialNum;
        
        % Loop through trials, starting from the second one (n-1 concept)
        for j = 2:height(session_data)
            current_trial = session_data(j, :);
            previous_trial = session_data(j-1, :);

            CurrentTrialNumber = current_trial.TrialNum;
            PrevTrialNumber = CurrentTrialNumber - 1;

            if ismember(PrevTrialNumber, TrialsIndices)
                Valids=Valids+1;
                % Check current and previous trial types
                if current_trial.EventValue == 4 && previous_trial.EventValue == 4
                    % TP and previous TP
                    TP_prev_TP = [TP_prev_TP; current_trial];
                    TP_prev_TP_Pre = [TP_prev_TP_Pre; previous_trial];
                elseif current_trial.EventValue == 4 && previous_trial.EventValue == 3
                    % TP and previous TA
                    TP_prev_TA = [TP_prev_TA; current_trial];
                    TP_prev_TA_Pre = [TP_prev_TA_Pre; previous_trial];
                elseif current_trial.EventValue == 3 && previous_trial.EventValue == 4
                    % TA and previous TP
                    TA_prev_TP = [TA_prev_TP; current_trial];
                    TA_prev_TP_Pre = [TA_prev_TP_Pre; previous_trial];
                elseif current_trial.EventValue == 3 && previous_trial.EventValue == 3
                    % TA and previous TA
                    TA_prev_TA = [TA_prev_TA; current_trial];
                    TA_prev_TA_Pre = [TA_prev_TA_Pre; previous_trial];
                end
            else
                counter = counter + 1;
            end
        end
        
        % Update the start index for the new session
        session_start_idx = i;
    end
    
    % Update prevTaskOrder
    prevTaskOrder = currentTaskOrder;
end

% Now TP_prev_TP, TP_prev_TA, TA_prev_TP, TA_prev_TA contain the trials
% from all sessions concatenated in the respective categories.



% Now TP_prev_TP, TP_prev_TA, TA_prev_TP, TA_prev_TA contain the trials
% from all sessions concatenated in the respective categories.
%% 

column_names = cell(1, 1600);

for i = 1:1600
    column_names{i} = ['bin', num2str(i)];
end



% TP-TP
TP_prev_TP_bins = TP_prev_TP{:, column_names(1:1600)};
TP_prev_TP_bins_mean=nanmean(TP_prev_TP_bins,1);
% TP-TA
TP_prev_TA_bins = TP_prev_TA{:, column_names(1:1600)};
TP_prev_TA_bins_mean=nanmean(TP_prev_TA_bins,1);
%TA-TP
TA_prev_TP_bins = TA_prev_TP{:, column_names(1:1600)};
TA_prev_TP_bins_mean=nanmean(TA_prev_TP_bins,1);
%TA-TA
TA_prev_TA_bins = TA_prev_TA{:, column_names(1:1600)};
TA_prev_TA_bins_mean=nanmean(TA_prev_TA_bins,1);


% Define time axis (assuming your data spans 1 second)
time_axis = linspace(-0.4, 1.2, 1600);

% Plot firing rates for each condition
figure;
hold on;
colors = get_distinguishable_colors(8); % You can also use other colormap functions like parula, jet, etc.
 
PlotPSTH(time_axis,TP_prev_TP_bins_mean,colors(1,:));
PlotPSTH(time_axis,TP_prev_TA_bins_mean,colors(2,:));

% Labeling
xlabel('Time (s)');
ylabel('Firing Rate (spikes/s)');
title('PSTH Firing Rates Across Conditions');
legend('TP Prev TP','TP prev TA');
ylim([40 90]);
grid on;
hold off;


figure;
hold on;
colors = get_distinguishable_colors(8); % You can also use other colormap functions like parula, jet, etc.
 
PlotPSTH(time_axis,TA_prev_TP_bins_mean,colors(1,:));
PlotPSTH(time_axis,TA_prev_TA_bins_mean,colors(2,:));

% Labeling
xlabel('Time (s)');
ylabel('Firing Rate (spikes/s)');
title('PSTH Firing Rates Across Conditions');
legend('TA Prev TP','TA prev TA');
ylim([40 90]);

grid on;
hold off;
%% 
TP_prev_TP_Final=TP_prev_TP;
TA_prev_TP_Final=TA_prev_TP;
TP_prev_TA_Final=TP_prev_TA;
TA_prev_TA_Final=TA_prev_TA;



% Get data for Monkey 1
Monkey3_TP_prev_TP_Final = TP_prev_TP_Final(TP_prev_TP_Final.Subject == 3, :);
Monkey3_TP_prev_TA_Final = TP_prev_TA_Final(TP_prev_TA_Final.Subject == 3, :);
Monkey3_TA_prev_TP_Final = TA_prev_TP_Final(TA_prev_TP_Final.Subject == 3, :);
Monkey3_TA_prev_TA_Final = TA_prev_TA_Final(TA_prev_TA_Final.Subject == 3, :);

% Get data for Monkey 4
Monkey4_TP_prev_TP_Final = TP_prev_TP_Final(TP_prev_TP_Final.Subject == 4, :);
Monkey4_TP_prev_TA_Final = TP_prev_TA_Final(TP_prev_TA_Final.Subject == 4, :);
Monkey4_TA_prev_TP_Final = TA_prev_TP_Final(TA_prev_TP_Final.Subject == 4, :);
Monkey4_TA_prev_TA_Final = TA_prev_TA_Final(TA_prev_TA_Final.Subject == 4, :);

% Define the columns (bins)
column_names = cell(1, 1600);
for i = 1:1600
    column_names{i} = ['bin', num2str(i)];
end

% Compute mean firing rates for Monkey 3
Monkey3_TP_prev_TP_Final_bins = Monkey3_TP_prev_TP_Final{:, column_names(1:1600)};
Monkey3_TP_prev_TP_Final_bins_mean = nanmean(Monkey3_TP_prev_TP_Final_bins, 1);

Monkey3_TP_prev_TA_Final_bins = Monkey3_TP_prev_TA_Final{:, column_names(1:1600)};
Monkey3_TP_prev_TA_Final_bins_mean = nanmean(Monkey3_TP_prev_TA_Final_bins, 1);

Monkey3_TA_prev_TP_Final_bins = Monkey3_TA_prev_TP_Final{:, column_names(1:1600)};
Monkey3_TA_prev_TP_Final_bins_mean = nanmean(Monkey3_TA_prev_TP_Final_bins, 1);

Monkey3_TA_prev_TA_Final_bins = Monkey3_TA_prev_TA_Final{:, column_names(1:1600)};
Monkey3_TA_prev_TA_Final_bins_mean = nanmean(Monkey3_TA_prev_TA_Final_bins, 1);

% Define the time axis
time_axis = linspace(-0.6, 1, 1600);

% Plot Monkey 1 firing rates
figure;
hold on;
colors = get_distinguishable_colors(8);
PlotPSTH(time_axis, Monkey3_TP_prev_TP_Final_bins_mean, colors(1,:));
PlotPSTH(time_axis, Monkey3_TP_prev_TA_Final_bins_mean, colors(2,:));

xlabel('Time (s)');
ylabel('Firing Rate (spikes/s)');
title('Monkey 3: PSTH Firing Rates Across Conditions');
legend('TP Prev TP', 'TP Prev TA');
grid on;
hold off;

figure;
hold on;
PlotPSTH(time_axis, Monkey3_TA_prev_TP_Final_bins_mean, colors(1,:));
PlotPSTH(time_axis, Monkey3_TA_prev_TA_Final_bins_mean, colors(2,:));

xlabel('Time (s)');
ylabel('Firing Rate (spikes/s)');
title('Monkey 3: PSTH Firing Rates Across Conditions');
legend('TA Prev TP', 'TA Prev TA');
grid on;
hold off;

% Repeat the same for Monkey 4
Monkey4_TP_prev_TP_Final_bins = Monkey4_TP_prev_TP_Final{:, column_names(1:1600)};
Monkey4_TP_prev_TP_Final_bins_mean = nanmean(Monkey4_TP_prev_TP_Final_bins, 1);

Monkey4_TP_prev_TA_Final_bins = Monkey4_TP_prev_TA_Final{:, column_names(1:1600)};
Monkey4_TP_prev_TA_Final_bins_mean = nanmean(Monkey4_TP_prev_TA_Final_bins, 1);

Monkey4_TA_prev_TP_Final_bins = Monkey4_TA_prev_TP_Final{:, column_names(1:1600)};
Monkey4_TA_prev_TP_Final_bins_mean = nanmean(Monkey4_TA_prev_TP_Final_bins, 1);

Monkey4_TA_prev_TA_Final_bins = Monkey4_TA_prev_TA_Final{:, column_names(1:1600)};
Monkey4_TA_prev_TA_Final_bins_mean = nanmean(Monkey4_TA_prev_TA_Final_bins, 1);

% Plot Monkey 4 firing rates
figure;
hold on;
PlotPSTH(time_axis, Monkey4_TP_prev_TP_Final_bins_mean, colors(1,:));
PlotPSTH(time_axis, Monkey4_TP_prev_TA_Final_bins_mean, colors(2,:));

xlabel('Time (s)');
ylabel('Firing Rate (spikes/s)');
title('Monkey 4: PSTH Firing Rates Across Conditions');
legend('TP Prev TP', 'TP Prev TA');
grid on;
hold off;

figure;
hold on;
PlotPSTH(time_axis, Monkey4_TA_prev_TP_Final_bins_mean, colors(1,:));
PlotPSTH(time_axis, Monkey4_TA_prev_TA_Final_bins_mean, colors(2,:));

xlabel('Time (s)');
ylabel('Firing Rate (spikes/s)');
title('Monkey 4: PSTH Firing Rates Across Conditions');
legend('TA Prev TP', 'TA Prev TA');
grid on;
hold off;






%% % Separate Eff and Ineff conditions 

% Separate Eff and Ineff conditions for TA_prev_TP
%TP - TP
Eff_TPTP = find(TP_prev_TP.Search_Type == 1);
Ineff_TPTP = find(TP_prev_TP.Search_Type == 0);

TP_prev_TP_Eff=TP_prev_TP(Eff_TPTP,:);
TP_prev_TP_Ineff=TP_prev_TP(Ineff_TPTP,:);

%TP - TA
Eff_TPTA = find(TP_prev_TA.Search_Type == 1);
Ineff_TPTA = find(TP_prev_TA.Search_Type == 0);

TP_prev_TA_Eff=TP_prev_TA(Eff_TPTA,:);
TP_prev_TA_Ineff=TP_prev_TA(Ineff_TPTA,:);

%TA - TP
Eff_TATP = find(TA_prev_TP.Search_Type == 1);
Ineff_TATP = find(TA_prev_TP.Search_Type == 0);

TA_prev_TP_Eff=TA_prev_TP(Eff_TATP,:);
TA_prev_TP_Ineff=TA_prev_TP(Ineff_TATP,:);

%TA - TA
Eff_TATA = find(TA_prev_TA.Search_Type == 1);
Ineff_TATA = find(TA_prev_TA.Search_Type == 0);

TA_prev_TA_Eff=TA_prev_TA(Eff_TATA,:);
TA_prev_TA_Ineff=TA_prev_TA(Ineff_TATA,:);


% Accessing to bins of PSTH
% TP-TP Eff
TP_prev_TP_Eff_bins = TP_prev_TP_Eff{:, column_names(1:1600)};
TP_prev_TP_Eff_bins_mean = nanmean(TP_prev_TP_Eff_bins, 1);

% TP-TP Ineff
TP_prev_TP_Ineff_bins = TP_prev_TP_Ineff{:, column_names(1:1600)};
TP_prev_TP_Ineff_bins_mean = nanmean(TP_prev_TP_Ineff_bins, 1);

% TP-TA Eff
TP_prev_TA_Eff_bins = TP_prev_TA_Eff{:, column_names(1:1600)};
TP_prev_TA_Eff_bins_mean = nanmean(TP_prev_TA_Eff_bins, 1);

% TP-TA Ineff
TP_prev_TA_Ineff_bins = TP_prev_TA_Ineff{:, column_names(1:1600)};
TP_prev_TA_Ineff_bins_mean = nanmean(TP_prev_TA_Ineff_bins, 1);

% TA-TP Eff
TA_prev_TP_Eff_bins = TA_prev_TP_Eff{:, column_names(1:1600)};
TA_prev_TP_Eff_bins_mean = nanmean(TA_prev_TP_Eff_bins, 1);

% TA-TP Ineff
TA_prev_TP_Ineff_bins = TA_prev_TP_Ineff{:, column_names(1:1600)};
TA_prev_TP_Ineff_bins_mean = nanmean(TA_prev_TP_Ineff_bins, 1);

% TA-TA Eff
TA_prev_TA_Eff_bins = TA_prev_TA_Eff{:, column_names(1:1600)};
TA_prev_TA_Eff_bins_mean = nanmean(TA_prev_TA_Eff_bins, 1);

% TA-TA Ineff
TA_prev_TA_Ineff_bins = TA_prev_TA_Ineff{:, column_names(1:1600)};
TA_prev_TA_Ineff_bins_mean = nanmean(TA_prev_TA_Ineff_bins, 1);





figure;
hold on;
colors = get_distinguishable_colors(8); % You can also use other colormap functions like parula, jet, etc.
 
PlotPSTH(time_axis,TP_prev_TP_Eff_bins_mean,colors(1,:));
PlotPSTH(time_axis,TP_prev_TA_Eff_bins_mean,colors(2,:));

% Labeling
xlabel('Time (s)');
ylabel('Firing Rate (spikes/s)');
title('PSTH Firing Rates Across Conditions for Effiecient');
legend('TP Prev TP','TP prev TA');
ylim([40 90]);
grid on;
hold off;


figure;
hold on;
colors = get_distinguishable_colors(8); % You can also use other colormap functions like parula, jet, etc.
 
PlotPSTH(time_axis,TP_prev_TP_Ineff_bins_mean,colors(1,:));
PlotPSTH(time_axis,TP_prev_TA_Ineff_bins_mean,colors(2,:));

% Labeling
xlabel('Time (s)');
ylabel('Firing Rate (spikes/s)');
title('PSTH Firing Rates Across Conditions for Inefficient');
legend('TP Prev TP','TP prev TA');
ylim([40 90]);

grid on;
hold off;

figure;
hold on;
colors = get_distinguishable_colors(8); % You can also use other colormap functions like parula, jet, etc.
 
PlotPSTH(time_axis,TA_prev_TP_Eff_bins_mean,colors(1,:));
PlotPSTH(time_axis,TA_prev_TA_Eff_bins_mean,colors(2,:));

% Labeling
xlabel('Time (s)');
ylabel('Firing Rate (spikes/s)');
title('PSTH Firing Rates Across Conditions Efficient');
legend('TA Prev TP','TA prev TA');
ylim([40 90]);
grid on;
hold off;


figure;
hold on;
colors = get_distinguishable_colors(8); % You can also use other colormap functions like parula, jet, etc.
 
PlotPSTH(time_axis,TA_prev_TP_Ineff_bins_mean,colors(1,:));
PlotPSTH(time_axis,TA_prev_TA_Ineff_bins_mean,colors(2,:));

% Labeling
xlabel('Time (s)');
ylabel('Firing Rate (spikes/s)');
title('PSTH Firing Rates Across Conditions for Inefficient');
legend('TA Prev TP','TA prev TA');
ylim([40 90]);

grid on;
hold off;


%% 
% Separate data by Subject (Monkey 3 and Monkey 4)
% Monkey 3
Monkey3_TP_prev_TP_Final = TP_prev_TP_Final(TP_prev_TP_Final.Subject == 3, :);
Monkey3_TP_prev_TA_Final = TP_prev_TA_Final(TP_prev_TA_Final.Subject == 3, :);
Monkey3_TA_prev_TP_Final = TA_prev_TP_Final(TA_prev_TP_Final.Subject == 3, :);
Monkey3_TA_prev_TA_Final = TA_prev_TA_Final(TA_prev_TA_Final.Subject == 3, :);

% Monkey 4
Monkey4_TP_prev_TP_Final = TP_prev_TP_Final(TP_prev_TP_Final.Subject == 4, :);
Monkey4_TP_prev_TA_Final = TP_prev_TA_Final(TP_prev_TA_Final.Subject == 4, :);
Monkey4_TA_prev_TP_Final = TA_prev_TP_Final(TA_prev_TP_Final.Subject == 4, :);
Monkey4_TA_prev_TA_Final = TA_prev_TA_Final(TA_prev_TA_Final.Subject == 4, :);

% For each monkey, separate Eff and Ineff conditions
% Monkey 3
Eff_TPTP_Monkey3 = find(Monkey3_TP_prev_TP_Final.Search_Type == 1);
Ineff_TPTP_Monkey3 = find(Monkey3_TP_prev_TP_Final.Search_Type == 0);
Monkey3_TP_prev_TP_Final_Eff = Monkey3_TP_prev_TP_Final(Eff_TPTP_Monkey3, :);
Monkey3_TP_prev_TP_Final_Ineff = Monkey3_TP_prev_TP_Final(Ineff_TPTP_Monkey3, :);

Eff_TPTA_Monkey3 = find(Monkey3_TP_prev_TA_Final.Search_Type == 1);
Ineff_TPTA_Monkey3 = find(Monkey3_TP_prev_TA_Final.Search_Type == 0);
Monkey3_TP_prev_TA_Final_Eff = Monkey3_TP_prev_TA_Final(Eff_TPTA_Monkey3, :);
Monkey3_TP_prev_TA_Final_Ineff = Monkey3_TP_prev_TA_Final(Ineff_TPTA_Monkey3, :);

Eff_TATP_Monkey3 = find(Monkey3_TA_prev_TP_Final.Search_Type == 1);
Ineff_TATP_Monkey3 = find(Monkey3_TA_prev_TP_Final.Search_Type == 0);
Monkey3_TA_prev_TP_Final_Eff = Monkey3_TA_prev_TP_Final(Eff_TATP_Monkey3, :);
Monkey3_TA_prev_TP_Final_Ineff = Monkey3_TA_prev_TP_Final(Ineff_TATP_Monkey3, :);

Eff_TATA_Monkey3 = find(Monkey3_TA_prev_TA_Final.Search_Type == 1);
Ineff_TATA_Monkey3 = find(Monkey3_TA_prev_TA_Final.Search_Type == 0);
Monkey3_TA_prev_TA_Final_Eff = Monkey3_TA_prev_TA_Final(Eff_TATA_Monkey3, :);
Monkey3_TA_prev_TA_Final_Ineff = Monkey3_TA_prev_TA_Final(Ineff_TATA_Monkey3, :);

% Monkey 4
Eff_TPTP_Monkey4 = find(Monkey4_TP_prev_TP_Final.Search_Type == 1);
Ineff_TPTP_Monkey4 = find(Monkey4_TP_prev_TP_Final.Search_Type == 0);
Monkey4_TP_prev_TP_Final_Eff = Monkey4_TP_prev_TP_Final(Eff_TPTP_Monkey4, :);
Monkey4_TP_prev_TP_Final_Ineff = Monkey4_TP_prev_TP_Final(Ineff_TPTP_Monkey4, :);

Eff_TPTA_Monkey4 = find(Monkey4_TP_prev_TA_Final.Search_Type == 1);
Ineff_TPTA_Monkey4 = find(Monkey4_TP_prev_TA_Final.Search_Type == 0);
Monkey4_TP_prev_TA_Final_Eff = Monkey4_TP_prev_TA_Final(Eff_TPTA_Monkey4, :);
Monkey4_TP_prev_TA_Final_Ineff = Monkey4_TP_prev_TA_Final(Ineff_TPTA_Monkey4, :);

Eff_TATP_Monkey4 = find(Monkey4_TA_prev_TP_Final.Search_Type == 1);
Ineff_TATP_Monkey4 = find(Monkey4_TA_prev_TP_Final.Search_Type == 0);
Monkey4_TA_prev_TP_Final_Eff = Monkey4_TA_prev_TP_Final(Eff_TATP_Monkey4, :);
Monkey4_TA_prev_TP_Final_Ineff = Monkey4_TA_prev_TP_Final(Ineff_TATP_Monkey4, :);

Eff_TATA_Monkey4 = find(Monkey4_TA_prev_TA_Final.Search_Type == 1);
Ineff_TATA_Monkey4 = find(Monkey4_TA_prev_TA_Final.Search_Type == 0);
Monkey4_TA_prev_TA_Final_Eff = Monkey4_TA_prev_TA_Final(Eff_TATA_Monkey4, :);
Monkey4_TA_prev_TA_Final_Ineff = Monkey4_TA_prev_TA_Final(Ineff_TATA_Monkey4, :);

% Now you can proceed with the plotting for each monkey and condition, similar to before.
% Accessing bins of PSTH for Monkey 3
% TP-TP Eff Monkey 3
Monkey3_TP_prev_TP_Final_Eff_bins = Monkey3_TP_prev_TP_Final_Eff{:, column_names(1:1600)};
Monkey3_TP_prev_TP_Final_Eff_bins_mean = nanmean(Monkey3_TP_prev_TP_Final_Eff_bins, 1);

% TP-TP Ineff Monkey 3
Monkey3_TP_prev_TP_Final_Ineff_bins = Monkey3_TP_prev_TP_Final_Ineff{:, column_names(1:1600)};
Monkey3_TP_prev_TP_Final_Ineff_bins_mean = nanmean(Monkey3_TP_prev_TP_Final_Ineff_bins, 1);

% TP-TA Eff Monkey 3
Monkey3_TP_prev_TA_Final_Eff_bins = Monkey3_TP_prev_TA_Final_Eff{:, column_names(1:1600)};
Monkey3_TP_prev_TA_Final_Eff_bins_mean = nanmean(Monkey3_TP_prev_TA_Final_Eff_bins, 1);

% TP-TA Ineff Monkey 3
Monkey3_TP_prev_TA_Final_Ineff_bins = Monkey3_TP_prev_TA_Final_Ineff{:, column_names(1:1600)};
Monkey3_TP_prev_TA_Final_Ineff_bins_mean = nanmean(Monkey3_TP_prev_TA_Final_Ineff_bins, 1);

% TA-TP Eff Monkey 3
Monkey3_TA_prev_TP_Final_Eff_bins = Monkey3_TA_prev_TP_Final_Eff{:, column_names(1:1600)};
Monkey3_TA_prev_TP_Final_Eff_bins_mean = nanmean(Monkey3_TA_prev_TP_Final_Eff_bins, 1);

% TA-TP Ineff Monkey 3
Monkey3_TA_prev_TP_Final_Ineff_bins = Monkey3_TA_prev_TP_Final_Ineff{:, column_names(1:1600)};
Monkey3_TA_prev_TP_Final_Ineff_bins_mean = nanmean(Monkey3_TA_prev_TP_Final_Ineff_bins, 1);

% TA-TA Eff Monkey 3
Monkey3_TA_prev_TA_Final_Eff_bins = Monkey3_TA_prev_TA_Final_Eff{:, column_names(1:1600)};
Monkey3_TA_prev_TA_Final_Eff_bins_mean = nanmean(Monkey3_TA_prev_TA_Final_Eff_bins, 1);

% TA-TA Ineff Monkey 3
Monkey3_TA_prev_TA_Final_Ineff_bins = Monkey3_TA_prev_TA_Final_Ineff{:, column_names(1:1600)};
Monkey3_TA_prev_TA_Final_Ineff_bins_mean = nanmean(Monkey3_TA_prev_TA_Final_Ineff_bins, 1);

% Now for Monkey 4
% TP-TP Eff Monkey 4
Monkey4_TP_prev_TP_Final_Eff_bins = Monkey4_TP_prev_TP_Final_Eff{:, column_names(1:1600)};
Monkey4_TP_prev_TP_Final_Eff_bins_mean = nanmean(Monkey4_TP_prev_TP_Final_Eff_bins, 1);

% TP-TP Ineff Monkey 4
Monkey4_TP_prev_TP_Final_Ineff_bins = Monkey4_TP_prev_TP_Final_Ineff{:, column_names(1:1600)};
Monkey4_TP_prev_TP_Final_Ineff_bins_mean = nanmean(Monkey4_TP_prev_TP_Final_Ineff_bins, 1);

% TP-TA Eff Monkey 4
Monkey4_TP_prev_TA_Final_Eff_bins = Monkey4_TP_prev_TA_Final_Eff{:, column_names(1:1600)};
Monkey4_TP_prev_TA_Final_Eff_bins_mean = nanmean(Monkey4_TP_prev_TA_Final_Eff_bins, 1);

% TP-TA Ineff Monkey 4
Monkey4_TP_prev_TA_Final_Ineff_bins = Monkey4_TP_prev_TA_Final_Ineff{:, column_names(1:1600)};
Monkey4_TP_prev_TA_Final_Ineff_bins_mean = nanmean(Monkey4_TP_prev_TA_Final_Ineff_bins, 1);

% TA-TP Eff Monkey 4
Monkey4_TA_prev_TP_Final_Eff_bins = Monkey4_TA_prev_TP_Final_Eff{:, column_names(1:1600)};
Monkey4_TA_prev_TP_Final_Eff_bins_mean = nanmean(Monkey4_TA_prev_TP_Final_Eff_bins, 1);

% TA-TP Ineff Monkey 4
Monkey4_TA_prev_TP_Final_Ineff_bins = Monkey4_TA_prev_TP_Final_Ineff{:, column_names(1:1600)};
Monkey4_TA_prev_TP_Final_Ineff_bins_mean = nanmean(Monkey4_TA_prev_TP_Final_Ineff_bins, 1);

% TA-TA Eff Monkey 4
Monkey4_TA_prev_TA_Final_Eff_bins = Monkey4_TA_prev_TA_Final_Eff{:, column_names(1:1600)};
Monkey4_TA_prev_TA_Final_Eff_bins_mean = nanmean(Monkey4_TA_prev_TA_Final_Eff_bins, 1);

% TA-TA Ineff Monkey 4
Monkey4_TA_prev_TA_Final_Ineff_bins = Monkey4_TA_prev_TA_Final_Ineff{:, column_names(1:1600)};
Monkey4_TA_prev_TA_Final_Ineff_bins_mean = nanmean(Monkey4_TA_prev_TA_Final_Ineff_bins, 1);

% Plot for Monkey 3
figure;
hold on;
colors = get_distinguishable_colors(8);
PlotPSTH(time_axis, Monkey3_TP_prev_TP_Final_Eff_bins_mean, colors(1,:));
PlotPSTH(time_axis, Monkey3_TP_prev_TA_Final_Eff_bins_mean, colors(2,:));
xlabel('Time (s)');
ylabel('Firing Rate (spikes/s)');
title('Monkey 3: PSTH Firing Rates Across Conditions for Efficient');
legend('TP Prev TP','TP Prev TA');
grid on;
hold off;

figure;
hold on;
PlotPSTH(time_axis, Monkey3_TP_prev_TP_Final_Ineff_bins_mean, colors(1,:));
PlotPSTH(time_axis, Monkey3_TP_prev_TA_Final_Ineff_bins_mean, colors(2,:));
xlabel('Time (s)');
ylabel('Firing Rate (spikes/s)');
title('Monkey 3: PSTH Firing Rates Across Conditions for Inefficient');
legend('TP Prev TP','TP Prev TA');
grid on;
hold off;

figure;
hold on;
PlotPSTH(time_axis, Monkey3_TA_prev_TP_Final_Eff_bins_mean, colors(1,:));
PlotPSTH(time_axis, Monkey3_TA_prev_TA_Final_Eff_bins_mean, colors(2,:));
xlabel('Time (s)');
ylabel('Firing Rate (spikes/s)');
title('Monkey 3: PSTH Firing Rates Across Conditions for Efficient');
legend('TA Prev TP','TA Prev TA');
grid on;
hold off;

figure;
hold on;
PlotPSTH(time_axis, Monkey3_TA_prev_TP_Final_Ineff_bins_mean, colors(1,:));
PlotPSTH(time_axis, Monkey3_TA_prev_TA_Final_Ineff_bins_mean, colors(2,:));
xlabel('Time (s)');
ylabel('Firing Rate (spikes/s)');
title('Monkey 3: PSTH Firing Rates Across Conditions for Inefficient');
legend('TA Prev TP','TA Prev TA');
grid on;
hold off;

% Plot for Monkey 4
figure;
hold on;
PlotPSTH(time_axis, Monkey4_TP_prev_TP_Final_Eff_bins_mean, colors(1,:));
PlotPSTH(time_axis, Monkey4_TP_prev_TA_Final_Eff_bins_mean, colors(2,:));
xlabel('Time (s)');
ylabel('Firing Rate (spikes/s)');
title('Monkey 4: PSTH Firing Rates Across Conditions for Efficient');
legend('TP Prev TP','TP Prev TA');
grid on;
hold off;

figure;
hold on;
PlotPSTH(time_axis, Monkey4_TP_prev_TP_Final_Ineff_bins_mean, colors(1,:));
PlotPSTH(time_axis, Monkey4_TP_prev_TA_Final_Ineff_bins_mean, colors(2,:));
xlabel('Time (s)');
ylabel('Firing Rate (spikes/s)');
title('Monkey 4: PSTH Firing Rates Across Conditions for Inefficient');
legend('TP Prev TP','TP Prev TA');
grid on;
hold off;

figure;
hold on;
PlotPSTH(time_axis, Monkey4_TA_prev_TP_Final_Eff_bins_mean, colors(1,:));
PlotPSTH(time_axis, Monkey4_TA_prev_TA_Final_Eff_bins_mean, colors(2,:));
xlabel('Time (s)');
ylabel('Firing Rate (spikes/s)');
title('Monkey 4: PSTH Firing Rates Across Conditions for Efficient');
legend('TA Prev TP','TA Prev TA');
grid on;
hold off;

figure;
hold on;
PlotPSTH(time_axis, Monkey4_TA_prev_TP_Final_Ineff_bins_mean, colors(1,:));
PlotPSTH(time_axis, Monkey4_TA_prev_TA_Final_Ineff_bins_mean, colors(2,:));
xlabel('Time (s)');
ylabel('Firing Rate (spikes/s)');
title('Monkey 4: PSTH Firing Rates Across Conditions for Inefficient');
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
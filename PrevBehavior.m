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
% C = Current P = Previous
% Conditions:
% C:TP/P:TP -- C:TP/P:TA -- C:TA/P:TP -- C:TA/P:TA


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



%% Current trial is TP

% TP_TP display size separation
DspSize3_TP = find(TP_prev_TP.DispSize == 3);
DspSize5_TP = find(TP_prev_TP.DispSize == 5);
DspSize7_TP = find(TP_prev_TP.DispSize == 7);
DspSize9_TP = find(TP_prev_TP.DispSize == 9);

TP_prev_TP_DspSize3 = TP_prev_TP(DspSize3_TP, :);
TP_prev_TP_DspSize5 = TP_prev_TP(DspSize5_TP, :);
TP_prev_TP_DspSize7 = TP_prev_TP(DspSize7_TP, :);
TP_prev_TP_DspSize9 = TP_prev_TP(DspSize9_TP, :);

% Search time Mean and SEM for TP_prev_TP
TP_prev_TP_DspSize3_SearchTime = nanmean(TP_prev_TP_DspSize3.SearchTime);
TP_prev_TP_DspSize5_SearchTime = nanmean(TP_prev_TP_DspSize5.SearchTime);
TP_prev_TP_DspSize7_SearchTime = nanmean(TP_prev_TP_DspSize7.SearchTime);
TP_prev_TP_DspSize9_SearchTime = nanmean(TP_prev_TP_DspSize9.SearchTime);

TP_prev_TP_DspSize3_SEM = nanstd(TP_prev_TP_DspSize3.SearchTime) / sqrt(length(TP_prev_TP_DspSize3.SearchTime));
TP_prev_TP_DspSize5_SEM = nanstd(TP_prev_TP_DspSize5.SearchTime) / sqrt(length(TP_prev_TP_DspSize5.SearchTime));
TP_prev_TP_DspSize7_SEM = nanstd(TP_prev_TP_DspSize7.SearchTime) / sqrt(length(TP_prev_TP_DspSize7.SearchTime));
TP_prev_TP_DspSize9_SEM = nanstd(TP_prev_TP_DspSize9.SearchTime) / sqrt(length(TP_prev_TP_DspSize9.SearchTime));

% TP_TA display size separation
DspSize3_TA = find(TP_prev_TA.DispSize == 3);
DspSize5_TA = find(TP_prev_TA.DispSize == 5);
DspSize7_TA = find(TP_prev_TA.DispSize == 7);
DspSize9_TA = find(TP_prev_TA.DispSize == 9);

TP_prev_TA_DspSize3 = TP_prev_TA(DspSize3_TA, :);
TP_prev_TA_DspSize5 = TP_prev_TA(DspSize5_TA, :);
TP_prev_TA_DspSize7 = TP_prev_TA(DspSize7_TA, :);
TP_prev_TA_DspSize9 = TP_prev_TA(DspSize9_TA, :);

% Search time Mean and SEM for TP_prev_TA
TP_prev_TA_DspSize3_SearchTime = nanmean(TP_prev_TA_DspSize3.SearchTime);
TP_prev_TA_DspSize5_SearchTime = nanmean(TP_prev_TA_DspSize5.SearchTime);
TP_prev_TA_DspSize7_SearchTime = nanmean(TP_prev_TA_DspSize7.SearchTime);
TP_prev_TA_DspSize9_SearchTime = nanmean(TP_prev_TA_DspSize9.SearchTime);

TP_prev_TA_DspSize3_SEM = nanstd(TP_prev_TA_DspSize3.SearchTime) / sqrt(length(TP_prev_TA_DspSize3.SearchTime));
TP_prev_TA_DspSize5_SEM = nanstd(TP_prev_TA_DspSize5.SearchTime) / sqrt(length(TP_prev_TA_DspSize5.SearchTime));
TP_prev_TA_DspSize7_SEM = nanstd(TP_prev_TA_DspSize7.SearchTime) / sqrt(length(TP_prev_TA_DspSize7.SearchTime));
TP_prev_TA_DspSize9_SEM = nanstd(TP_prev_TA_DspSize9.SearchTime) / sqrt(length(TP_prev_TA_DspSize9.SearchTime));

% Prepare data for plotting
display_sizes = [3, 5, 7, 9];

% Mean search times for TP_prev_TP and TP_prev_TA
mean_search_times_TP = [TP_prev_TP_DspSize3_SearchTime, TP_prev_TP_DspSize5_SearchTime, ...
                        TP_prev_TP_DspSize7_SearchTime, TP_prev_TP_DspSize9_SearchTime];
mean_search_times_TA = [TP_prev_TA_DspSize3_SearchTime, TP_prev_TA_DspSize5_SearchTime, ...
                        TP_prev_TA_DspSize7_SearchTime, TP_prev_TA_DspSize9_SearchTime];

% SEM search times for TP_prev_TP and TP_prev_TA
sem_search_times_TP = [TP_prev_TP_DspSize3_SEM, TP_prev_TP_DspSize5_SEM, ...
                       TP_prev_TP_DspSize7_SEM, TP_prev_TP_DspSize9_SEM];
sem_search_times_TA = [TP_prev_TA_DspSize3_SEM, TP_prev_TA_DspSize5_SEM, ...
                       TP_prev_TA_DspSize7_SEM, TP_prev_TA_DspSize9_SEM];

% Plot with error bars (mean + SEM) for both states
figure;
errorbar(display_sizes, mean_search_times_TP, sem_search_times_TP, 'o-', 'LineWidth', 2, 'DisplayName', 'TP Prev TP');
hold on;
errorbar(display_sizes, mean_search_times_TA, sem_search_times_TA, 's-', 'LineWidth', 2, 'DisplayName', 'TP Prev TA');
set(gca,'FontSize',23); 

% Add labels and title
xlabel('Display Size');
ylabel('Mean Search Time');
title('Mean Search Time by Display Size for TP Prev TP and TP Prev TA');
legend('show');
xlim([2,10]);
grid on;
hold off;




%% Monkey Separate TP current trial
% Current trial is TP with Subject-based separation (Subjects 3 and 4)

% Separate data by Subject (Monkey 3 and Monkey 4)
% Monkey 3
Monkey3_TP_prev_TP_Final = TP_prev_TP_Final(TP_prev_TP_Final.Subject == 3, :);
Monkey3_TP_prev_TA_Final = TP_prev_TA_Final(TP_prev_TA_Final.Subject == 3, :);

% Monkey 4
Monkey4_TP_prev_TP_Final = TP_prev_TP_Final(TP_prev_TP_Final.Subject == 4, :);
Monkey4_TP_prev_TA_Final = TP_prev_TA_Final(TP_prev_TA_Final.Subject == 4, :);

% For Monkey 3
% TP_TP display size separation for Monkey 3
DspSize3_TP_Monkey3 = find(Monkey3_TP_prev_TP_Final.DispSize == 3);
DspSize5_TP_Monkey3 = find(Monkey3_TP_prev_TP_Final.DispSize == 5);
DspSize7_TP_Monkey3 = find(Monkey3_TP_prev_TP_Final.DispSize == 7);
DspSize9_TP_Monkey3 = find(Monkey3_TP_prev_TP_Final.DispSize == 9);

Monkey3_TP_prev_TP_Final_DspSize3 = Monkey3_TP_prev_TP_Final(DspSize3_TP_Monkey3, :);
Monkey3_TP_prev_TP_Final_DspSize5 = Monkey3_TP_prev_TP_Final(DspSize5_TP_Monkey3, :);
Monkey3_TP_prev_TP_Final_DspSize7 = Monkey3_TP_prev_TP_Final(DspSize7_TP_Monkey3, :);
Monkey3_TP_prev_TP_Final_DspSize9 = Monkey3_TP_prev_TP_Final(DspSize9_TP_Monkey3, :);

% Search time Mean and SEM for TP_prev_TP_Final for Monkey 3
Monkey3_TP_prev_TP_Final_DspSize3_SearchTime = nanmean(Monkey3_TP_prev_TP_Final_DspSize3.SearchTime);
Monkey3_TP_prev_TP_Final_DspSize5_SearchTime = nanmean(Monkey3_TP_prev_TP_Final_DspSize5.SearchTime);
Monkey3_TP_prev_TP_Final_DspSize7_SearchTime = nanmean(Monkey3_TP_prev_TP_Final_DspSize7.SearchTime);
Monkey3_TP_prev_TP_Final_DspSize9_SearchTime = nanmean(Monkey3_TP_prev_TP_Final_DspSize9.SearchTime);

Monkey3_TP_prev_TP_Final_DspSize3_SEM = nanstd(Monkey3_TP_prev_TP_Final_DspSize3.SearchTime) / sqrt(length(Monkey3_TP_prev_TP_Final_DspSize3.SearchTime));
Monkey3_TP_prev_TP_Final_DspSize5_SEM = nanstd(Monkey3_TP_prev_TP_Final_DspSize5.SearchTime) / sqrt(length(Monkey3_TP_prev_TP_Final_DspSize5.SearchTime));
Monkey3_TP_prev_TP_Final_DspSize7_SEM = nanstd(Monkey3_TP_prev_TP_Final_DspSize7.SearchTime) / sqrt(length(Monkey3_TP_prev_TP_Final_DspSize7.SearchTime));
Monkey3_TP_prev_TP_Final_DspSize9_SEM = nanstd(Monkey3_TP_prev_TP_Final_DspSize9.SearchTime) / sqrt(length(Monkey3_TP_prev_TP_Final_DspSize9.SearchTime));

% Repeat the same for TP_TA display size separation for Monkey 3
DspSize3_TA_Monkey3 = find(Monkey3_TP_prev_TA_Final.DispSize == 3);
DspSize5_TA_Monkey3 = find(Monkey3_TP_prev_TA_Final.DispSize == 5);
DspSize7_TA_Monkey3 = find(Monkey3_TP_prev_TA_Final.DispSize == 7);
DspSize9_TA_Monkey3 = find(Monkey3_TP_prev_TA_Final.DispSize == 9);

Monkey3_TP_prev_TA_Final_DspSize3 = Monkey3_TP_prev_TA_Final(DspSize3_TA_Monkey3, :);
Monkey3_TP_prev_TA_Final_DspSize5 = Monkey3_TP_prev_TA_Final(DspSize5_TA_Monkey3, :);
Monkey3_TP_prev_TA_Final_DspSize7 = Monkey3_TP_prev_TA_Final(DspSize7_TA_Monkey3, :);
Monkey3_TP_prev_TA_Final_DspSize9 = Monkey3_TP_prev_TA_Final(DspSize9_TA_Monkey3, :);

% Search time Mean and SEM for TP_prev_TA_Final for Monkey 3
Monkey3_TP_prev_TA_Final_DspSize3_SearchTime = nanmean(Monkey3_TP_prev_TA_Final_DspSize3.SearchTime);
Monkey3_TP_prev_TA_Final_DspSize5_SearchTime = nanmean(Monkey3_TP_prev_TA_Final_DspSize5.SearchTime);
Monkey3_TP_prev_TA_Final_DspSize7_SearchTime = nanmean(Monkey3_TP_prev_TA_Final_DspSize7.SearchTime);
Monkey3_TP_prev_TA_Final_DspSize9_SearchTime = nanmean(Monkey3_TP_prev_TA_Final_DspSize9.SearchTime);

Monkey3_TP_prev_TA_Final_DspSize3_SEM = nanstd(Monkey3_TP_prev_TA_Final_DspSize3.SearchTime) / sqrt(length(Monkey3_TP_prev_TA_Final_DspSize3.SearchTime));
Monkey3_TP_prev_TA_Final_DspSize5_SEM = nanstd(Monkey3_TP_prev_TA_Final_DspSize5.SearchTime) / sqrt(length(Monkey3_TP_prev_TA_Final_DspSize5.SearchTime));
Monkey3_TP_prev_TA_Final_DspSize7_SEM = nanstd(Monkey3_TP_prev_TA_Final_DspSize7.SearchTime) / sqrt(length(Monkey3_TP_prev_TA_Final_DspSize7.SearchTime));
Monkey3_TP_prev_TA_Final_DspSize9_SEM = nanstd(Monkey3_TP_prev_TA_Final_DspSize9.SearchTime) / sqrt(length(Monkey3_TP_prev_TA_Final_DspSize9.SearchTime));

% For Monkey 4
% TP_TP display size separation for Monkey 4
DspSize3_TP_Monkey4 = find(Monkey4_TP_prev_TP_Final.DispSize == 3);
DspSize5_TP_Monkey4 = find(Monkey4_TP_prev_TP_Final.DispSize == 5);
DspSize7_TP_Monkey4 = find(Monkey4_TP_prev_TP_Final.DispSize == 7);
DspSize9_TP_Monkey4 = find(Monkey4_TP_prev_TP_Final.DispSize == 9);

Monkey4_TP_prev_TP_Final_DspSize3 = Monkey4_TP_prev_TP_Final(DspSize3_TP_Monkey4, :);
Monkey4_TP_prev_TP_Final_DspSize5 = Monkey4_TP_prev_TP_Final(DspSize5_TP_Monkey4, :);
Monkey4_TP_prev_TP_Final_DspSize7 = Monkey4_TP_prev_TP_Final(DspSize7_TP_Monkey4, :);
Monkey4_TP_prev_TP_Final_DspSize9 = Monkey4_TP_prev_TP_Final(DspSize9_TP_Monkey4, :);

% Search time Mean and SEM for TP_prev_TP_Final for Monkey 4
Monkey4_TP_prev_TP_Final_DspSize3_SearchTime = nanmean(Monkey4_TP_prev_TP_Final_DspSize3.SearchTime);
Monkey4_TP_prev_TP_Final_DspSize5_SearchTime = nanmean(Monkey4_TP_prev_TP_Final_DspSize5.SearchTime);
Monkey4_TP_prev_TP_Final_DspSize7_SearchTime = nanmean(Monkey4_TP_prev_TP_Final_DspSize7.SearchTime);
Monkey4_TP_prev_TP_Final_DspSize9_SearchTime = nanmean(Monkey4_TP_prev_TP_Final_DspSize9.SearchTime);

Monkey4_TP_prev_TP_Final_DspSize3_SEM = nanstd(Monkey4_TP_prev_TP_Final_DspSize3.SearchTime) / sqrt(length(Monkey4_TP_prev_TP_Final_DspSize3.SearchTime));
Monkey4_TP_prev_TP_Final_DspSize5_SEM = nanstd(Monkey4_TP_prev_TP_Final_DspSize5.SearchTime) / sqrt(length(Monkey4_TP_prev_TP_Final_DspSize5.SearchTime));
Monkey4_TP_prev_TP_Final_DspSize7_SEM = nanstd(Monkey4_TP_prev_TP_Final_DspSize7.SearchTime) / sqrt(length(Monkey4_TP_prev_TP_Final_DspSize7.SearchTime));
Monkey4_TP_prev_TP_Final_DspSize9_SEM = nanstd(Monkey4_TP_prev_TP_Final_DspSize9.SearchTime) / sqrt(length(Monkey4_TP_prev_TP_Final_DspSize9.SearchTime));

% Repeat the same for TP_TA display size separation for Monkey 4
DspSize3_TA_Monkey4 = find(Monkey4_TP_prev_TA_Final.DispSize == 3);
DspSize5_TA_Monkey4 = find(Monkey4_TP_prev_TA_Final.DispSize == 5);
DspSize7_TA_Monkey4 = find(Monkey4_TP_prev_TA_Final.DispSize == 7);
DspSize9_TA_Monkey4 = find(Monkey4_TP_prev_TA_Final.DispSize == 9);

Monkey4_TP_prev_TA_Final_DspSize3 = Monkey4_TP_prev_TA_Final(DspSize3_TA_Monkey4, :);
Monkey4_TP_prev_TA_Final_DspSize5 = Monkey4_TP_prev_TA_Final(DspSize5_TA_Monkey4, :);
Monkey4_TP_prev_TA_Final_DspSize7 = Monkey4_TP_prev_TA_Final(DspSize7_TA_Monkey4, :);
Monkey4_TP_prev_TA_Final_DspSize9 = Monkey4_TP_prev_TA_Final(DspSize9_TA_Monkey4, :);

% Search time Mean and SEM for TP_prev_TA_Final for Monkey 4
Monkey4_TP_prev_TA_Final_DspSize3_SearchTime = nanmean(Monkey4_TP_prev_TA_Final_DspSize3.SearchTime);
Monkey4_TP_prev_TA_Final_DspSize5_SearchTime = nanmean(Monkey4_TP_prev_TA_Final_DspSize5.SearchTime);
Monkey4_TP_prev_TA_Final_DspSize7_SearchTime = nanmean(Monkey4_TP_prev_TA_Final_DspSize7.SearchTime);
Monkey4_TP_prev_TA_Final_DspSize9_SearchTime = nanmean(Monkey4_TP_prev_TA_Final_DspSize9.SearchTime);

Monkey4_TP_prev_TA_Final_DspSize3_SEM = nanstd(Monkey4_TP_prev_TA_Final_DspSize3.SearchTime) / sqrt(length(Monkey4_TP_prev_TA_Final_DspSize3.SearchTime));
Monkey4_TP_prev_TA_Final_DspSize5_SEM = nanstd(Monkey4_TP_prev_TA_Final_DspSize5.SearchTime) / sqrt(length(Monkey4_TP_prev_TA_Final_DspSize5.SearchTime));
Monkey4_TP_prev_TA_Final_DspSize7_SEM = nanstd(Monkey4_TP_prev_TA_Final_DspSize7.SearchTime) / sqrt(length(Monkey4_TP_prev_TA_Final_DspSize7.SearchTime));
Monkey4_TP_prev_TA_Final_DspSize9_SEM = nanstd(Monkey4_TP_prev_TA_Final_DspSize9.SearchTime) / sqrt(length(Monkey4_TP_prev_TA_Final_DspSize9.SearchTime));

% Prepare data for plotting
display_sizes = [3, 5, 7, 9];

% Mean search times for Monkey 3 and Monkey 4
mean_search_times_TP_Monkey3 = [Monkey3_TP_prev_TP_Final_DspSize3_SearchTime, Monkey3_TP_prev_TP_Final_DspSize5_SearchTime, ...
                                Monkey3_TP_prev_TP_Final_DspSize7_SearchTime, Monkey3_TP_prev_TP_Final_DspSize9_SearchTime];
mean_search_times_TA_Monkey3 = [Monkey3_TP_prev_TA_Final_DspSize3_SearchTime, Monkey3_TP_prev_TA_Final_DspSize5_SearchTime, ...
                                Monkey3_TP_prev_TA_Final_DspSize7_SearchTime, Monkey3_TP_prev_TA_Final_DspSize9_SearchTime];

mean_search_times_TP_Monkey4 = [Monkey4_TP_prev_TP_Final_DspSize3_SearchTime, Monkey4_TP_prev_TP_Final_DspSize5_SearchTime, ...
                                Monkey4_TP_prev_TP_Final_DspSize7_SearchTime, Monkey4_TP_prev_TP_Final_DspSize9_SearchTime];
mean_search_times_TA_Monkey4 = [Monkey4_TP_prev_TA_Final_DspSize3_SearchTime, Monkey4_TP_prev_TA_Final_DspSize5_SearchTime, ...
                                Monkey4_TP_prev_TA_Final_DspSize7_SearchTime, Monkey4_TP_prev_TA_Final_DspSize9_SearchTime];

% SEM search times for Monkey 3 and Monkey 4
sem_search_times_TP_Monkey3 = [Monkey3_TP_prev_TP_Final_DspSize3_SEM, Monkey3_TP_prev_TP_Final_DspSize5_SEM, ...
                               Monkey3_TP_prev_TP_Final_DspSize7_SEM, Monkey3_TP_prev_TP_Final_DspSize9_SEM];
sem_search_times_TA_Monkey3 = [Monkey3_TP_prev_TA_Final_DspSize3_SEM, Monkey3_TP_prev_TA_Final_DspSize5_SEM, ...
                               Monkey3_TP_prev_TA_Final_DspSize7_SEM, Monkey3_TP_prev_TA_Final_DspSize9_SEM];

sem_search_times_TP_Monkey4 = [Monkey4_TP_prev_TP_Final_DspSize3_SEM, Monkey4_TP_prev_TP_Final_DspSize5_SEM, ...
                               Monkey4_TP_prev_TP_Final_DspSize7_SEM, Monkey4_TP_prev_TP_Final_DspSize9_SEM];
sem_search_times_TA_Monkey4 = [Monkey4_TP_prev_TA_Final_DspSize3_SEM, Monkey4_TP_prev_TA_Final_DspSize5_SEM, ...
                               Monkey4_TP_prev_TA_Final_DspSize7_SEM, Monkey4_TP_prev_TA_Final_DspSize9_SEM];

% Plot with error bars for both Monkeys and conditions

% Monkey 3
figure;
errorbar(display_sizes, mean_search_times_TP_Monkey3, sem_search_times_TP_Monkey3, 's-', 'LineWidth', 2, 'DisplayName', 'Monkey 3: TP Prev TP', 'Color', 'r');
hold on;
errorbar(display_sizes, mean_search_times_TA_Monkey3, sem_search_times_TA_Monkey3, 'o-', 'LineWidth', 2, 'DisplayName', 'Monkey 3: TP Prev TA', 'Color', 'b');
xlabel('Display Size');
ylabel('Mean Search Time');
title('Mean Search Time by Display Size for Monkey 3 (TP)');
legend('show');
grid on;
hold off;

% Monkey 4
figure;
errorbar(display_sizes, mean_search_times_TP_Monkey4, sem_search_times_TP_Monkey4, 's-', 'LineWidth', 2, 'DisplayName', 'Monkey 4: TP Prev TP', 'Color', 'r');
hold on;
errorbar(display_sizes, mean_search_times_TA_Monkey4, sem_search_times_TA_Monkey4, 'o-', 'LineWidth', 2, 'DisplayName', 'Monkey 4: TP Prev TA', 'Color', 'b');
xlabel('Display Size');
ylabel('Mean Search Time');
title('Mean Search Time by Display Size for Monkey 4 (TP)');
legend('show');
grid on;
hold off;






%% Current trial is TA
% TA_TP display size separation
DspSize3_TP = find(TA_prev_TP.DispSize == 3);
DspSize5_TP = find(TA_prev_TP.DispSize == 5);
DspSize7_TP = find(TA_prev_TP.DispSize == 7);
DspSize9_TP = find(TA_prev_TP.DispSize == 9);

TA_prev_TP_DspSize3 = TA_prev_TP(DspSize3_TP, :);
TA_prev_TP_DspSize5 = TA_prev_TP(DspSize5_TP, :);
TA_prev_TP_DspSize7 = TA_prev_TP(DspSize7_TP, :);
TA_prev_TP_DspSize9 = TA_prev_TP(DspSize9_TP, :);

% Search time Mean and SEM for TA_prev_TP
TA_prev_TP_DspSize3_SearchTime = nanmean(TA_prev_TP_DspSize3.SearchTime);
TA_prev_TP_DspSize5_SearchTime = nanmean(TA_prev_TP_DspSize5.SearchTime);
TA_prev_TP_DspSize7_SearchTime = nanmean(TA_prev_TP_DspSize7.SearchTime);
TA_prev_TP_DspSize9_SearchTime = nanmean(TA_prev_TP_DspSize9.SearchTime);

TA_prev_TP_DspSize3_SEM = nanstd(TA_prev_TP_DspSize3.SearchTime) / sqrt(length(TA_prev_TP_DspSize3.SearchTime));
TA_prev_TP_DspSize5_SEM = nanstd(TA_prev_TP_DspSize5.SearchTime) / sqrt(length(TA_prev_TP_DspSize5.SearchTime));
TA_prev_TP_DspSize7_SEM = nanstd(TA_prev_TP_DspSize7.SearchTime) / sqrt(length(TA_prev_TP_DspSize7.SearchTime));
TA_prev_TP_DspSize9_SEM = nanstd(TA_prev_TP_DspSize9.SearchTime) / sqrt(length(TA_prev_TP_DspSize9.SearchTime));

% TA_TA display size separation
DspSize3_TA = find(TA_prev_TA.DispSize == 3);
DspSize5_TA = find(TA_prev_TA.DispSize == 5);
DspSize7_TA = find(TA_prev_TA.DispSize == 7);
DspSize9_TA = find(TA_prev_TA.DispSize == 9);

TA_prev_TA_DspSize3 = TA_prev_TA(DspSize3_TA, :);
TA_prev_TA_DspSize5 = TA_prev_TA(DspSize5_TA, :);
TA_prev_TA_DspSize7 = TA_prev_TA(DspSize7_TA, :);
TA_prev_TA_DspSize9 = TA_prev_TA(DspSize9_TA, :);

% Search time Mean and SEM for TA_prev_TA
TA_prev_TA_DspSize3_SearchTime = nanmean(TA_prev_TA_DspSize3.SearchTime);
TA_prev_TA_DspSize5_SearchTime = nanmean(TA_prev_TA_DspSize5.SearchTime);
TA_prev_TA_DspSize7_SearchTime = nanmean(TA_prev_TA_DspSize7.SearchTime);
TA_prev_TA_DspSize9_SearchTime = nanmean(TA_prev_TA_DspSize9.SearchTime);

TA_prev_TA_DspSize3_SEM = nanstd(TA_prev_TA_DspSize3.SearchTime) / sqrt(length(TA_prev_TA_DspSize3.SearchTime));
TA_prev_TA_DspSize5_SEM = nanstd(TA_prev_TA_DspSize5.SearchTime) / sqrt(length(TA_prev_TA_DspSize5.SearchTime));
TA_prev_TA_DspSize7_SEM = nanstd(TA_prev_TA_DspSize7.SearchTime) / sqrt(length(TA_prev_TA_DspSize7.SearchTime));
TA_prev_TA_DspSize9_SEM = nanstd(TA_prev_TA_DspSize9.SearchTime) / sqrt(length(TA_prev_TA_DspSize9.SearchTime));

% Prepare data for plotting
display_sizes = [3, 5, 7, 9];

% Mean search times for TA_prev_TP and TA_prev_TA
mean_search_times_TP = [TA_prev_TP_DspSize3_SearchTime, TA_prev_TP_DspSize5_SearchTime, ...
                        TA_prev_TP_DspSize7_SearchTime, TA_prev_TP_DspSize9_SearchTime];
mean_search_times_TA = [TA_prev_TA_DspSize3_SearchTime, TA_prev_TA_DspSize5_SearchTime, ...
                        TA_prev_TA_DspSize7_SearchTime, TA_prev_TA_DspSize9_SearchTime];

% SEM search times for TA_prev_TP and TA_prev_TA
sem_search_times_TP = [TA_prev_TP_DspSize3_SEM, TA_prev_TP_DspSize5_SEM, ...
                       TA_prev_TP_DspSize7_SEM, TA_prev_TP_DspSize9_SEM];
sem_search_times_TA = [TA_prev_TA_DspSize3_SEM, TA_prev_TA_DspSize5_SEM, ...
                       TA_prev_TA_DspSize7_SEM, TA_prev_TA_DspSize9_SEM];

% Plot with error bars (mean + SEM) for both states
figure;
errorbar(display_sizes, mean_search_times_TP, sem_search_times_TP, 'o-', 'LineWidth', 2, 'DisplayName', 'TA Prev TP');
hold on;
errorbar(display_sizes, mean_search_times_TA, sem_search_times_TA, 's-', 'LineWidth', 2, 'DisplayName', 'TA Prev TA');
set(gca,'FontSize',23); 

% Add labels and title
xlabel('Display Size');
ylabel('Mean Search Time');
title('Mean Search Time by Display Size for TA Prev: TP and TA Prev: TA');
legend('show');
xlim([2,10]);

grid on;
hold off;




%% Current trial TP - EFF/INEFF separated
% Separate Eff and Ineff conditions for TP_prev_TP
Eff_TP = find(TP_prev_TP.Search_Type == 1);
Ineff_TP = find(TP_prev_TP.Search_Type == 0);

TP_prev_TP_Eff = TP_prev_TP(Eff_TP, :);
TP_prev_TP_Ineff = TP_prev_TP(Ineff_TP, :);


% TP_TP display size separation for Eff
DspSize3_TP_Eff = find(TP_prev_TP_Eff.DispSize == 3);
DspSize5_TP_Eff = find(TP_prev_TP_Eff.DispSize == 5);
DspSize7_TP_Eff = find(TP_prev_TP_Eff.DispSize == 7);
DspSize9_TP_Eff = find(TP_prev_TP_Eff.DispSize == 9);

% Mean and SEM for Eff condition
TP_prev_TP_Eff_DspSize3_SearchTime = nanmean(TP_prev_TP_Eff(DspSize3_TP_Eff, :).SearchTime);
TP_prev_TP_Eff_DspSize5_SearchTime = nanmean(TP_prev_TP_Eff(DspSize5_TP_Eff, :).SearchTime);
TP_prev_TP_Eff_DspSize7_SearchTime = nanmean(TP_prev_TP_Eff(DspSize7_TP_Eff, :).SearchTime);
TP_prev_TP_Eff_DspSize9_SearchTime = nanmean(TP_prev_TP_Eff(DspSize9_TP_Eff, :).SearchTime);

TP_prev_TP_Eff_DspSize3_SEM = nanstd(TP_prev_TP_Eff(DspSize3_TP_Eff, :).SearchTime) / sqrt(length(DspSize3_TP_Eff));
TP_prev_TP_Eff_DspSize5_SEM = nanstd(TP_prev_TP_Eff(DspSize5_TP_Eff, :).SearchTime) / sqrt(length(DspSize5_TP_Eff));
TP_prev_TP_Eff_DspSize7_SEM = nanstd(TP_prev_TP_Eff(DspSize7_TP_Eff, :).SearchTime) / sqrt(length(DspSize7_TP_Eff));
TP_prev_TP_Eff_DspSize9_SEM = nanstd(TP_prev_TP_Eff(DspSize9_TP_Eff, :).SearchTime) / sqrt(length(DspSize9_TP_Eff));

% TP_TP display size separation for Ineff
DspSize3_TP_Ineff = find(TP_prev_TP_Ineff.DispSize == 3);
DspSize5_TP_Ineff = find(TP_prev_TP_Ineff.DispSize == 5);
DspSize7_TP_Ineff = find(TP_prev_TP_Ineff.DispSize == 7);
DspSize9_TP_Ineff = find(TP_prev_TP_Ineff.DispSize == 9);

% Mean and SEM for Ineff condition
TP_prev_TP_Ineff_DspSize3_SearchTime = nanmean(TP_prev_TP_Ineff(DspSize3_TP_Ineff, :).SearchTime);
TP_prev_TP_Ineff_DspSize5_SearchTime = nanmean(TP_prev_TP_Ineff(DspSize5_TP_Ineff, :).SearchTime);
TP_prev_TP_Ineff_DspSize7_SearchTime = nanmean(TP_prev_TP_Ineff(DspSize7_TP_Ineff, :).SearchTime);
TP_prev_TP_Ineff_DspSize9_SearchTime = nanmean(TP_prev_TP_Ineff(DspSize9_TP_Ineff, :).SearchTime);

TP_prev_TP_Ineff_DspSize3_SEM = nanstd(TP_prev_TP_Ineff(DspSize3_TP_Ineff, :).SearchTime) / sqrt(length(DspSize3_TP_Ineff));
TP_prev_TP_Ineff_DspSize5_SEM = nanstd(TP_prev_TP_Ineff(DspSize5_TP_Ineff, :).SearchTime) / sqrt(length(DspSize5_TP_Ineff));
TP_prev_TP_Ineff_DspSize7_SEM = nanstd(TP_prev_TP_Ineff(DspSize7_TP_Ineff, :).SearchTime) / sqrt(length(DspSize7_TP_Ineff));
TP_prev_TP_Ineff_DspSize9_SEM = nanstd(TP_prev_TP_Ineff(DspSize9_TP_Ineff, :).SearchTime) / sqrt(length(DspSize9_TP_Ineff));


% Separate Eff and Ineff conditions for TP_prev_TA
Eff_TA = find(TP_prev_TA.Search_Type == 1);
Ineff_TA = find(TP_prev_TA.Search_Type == 0);

TP_prev_TA_Eff = TP_prev_TA(Eff_TA, :);
TP_prev_TA_Ineff = TP_prev_TA(Ineff_TA, :);

% TP_TA display size separation for Eff
DspSize3_TA_Eff = find(TP_prev_TA_Eff.DispSize == 3);
DspSize5_TA_Eff = find(TP_prev_TA_Eff.DispSize == 5);
DspSize7_TA_Eff = find(TP_prev_TA_Eff.DispSize == 7);
DspSize9_TA_Eff = find(TP_prev_TA_Eff.DispSize == 9);

% Mean and SEM for Eff condition
TP_prev_TA_Eff_DspSize3_SearchTime = nanmean(TP_prev_TA_Eff(DspSize3_TA_Eff, :).SearchTime);
TP_prev_TA_Eff_DspSize5_SearchTime = nanmean(TP_prev_TA_Eff(DspSize5_TA_Eff, :).SearchTime);
TP_prev_TA_Eff_DspSize7_SearchTime = nanmean(TP_prev_TA_Eff(DspSize7_TA_Eff, :).SearchTime);
TP_prev_TA_Eff_DspSize9_SearchTime = nanmean(TP_prev_TA_Eff(DspSize9_TA_Eff, :).SearchTime);

TP_prev_TA_Eff_DspSize3_SEM = nanstd(TP_prev_TA_Eff(DspSize3_TA_Eff, :).SearchTime) / sqrt(length(DspSize3_TA_Eff));
TP_prev_TA_Eff_DspSize5_SEM = nanstd(TP_prev_TA_Eff(DspSize5_TA_Eff, :).SearchTime) / sqrt(length(DspSize5_TA_Eff));
TP_prev_TA_Eff_DspSize7_SEM = nanstd(TP_prev_TA_Eff(DspSize7_TA_Eff, :).SearchTime) / sqrt(length(DspSize7_TA_Eff));
TP_prev_TA_Eff_DspSize9_SEM = nanstd(TP_prev_TA_Eff(DspSize9_TA_Eff, :).SearchTime) / sqrt(length(DspSize9_TA_Eff));

% TP_TA display size separation for Ineff
DspSize3_TA_Ineff = find(TP_prev_TA_Ineff.DispSize == 3);
DspSize5_TA_Ineff = find(TP_prev_TA_Ineff.DispSize == 5);
DspSize7_TA_Ineff = find(TP_prev_TA_Ineff.DispSize == 7);
DspSize9_TA_Ineff = find(TP_prev_TA_Ineff.DispSize == 9);

% Mean and SEM for Ineff condition
TP_prev_TA_Ineff_DspSize3_SearchTime = nanmean(TP_prev_TA_Ineff(DspSize3_TA_Ineff, :).SearchTime);
TP_prev_TA_Ineff_DspSize5_SearchTime = nanmean(TP_prev_TA_Ineff(DspSize5_TA_Ineff, :).SearchTime);
TP_prev_TA_Ineff_DspSize7_SearchTime = nanmean(TP_prev_TA_Ineff(DspSize7_TA_Ineff, :).SearchTime);
TP_prev_TA_Ineff_DspSize9_SearchTime = nanmean(TP_prev_TA_Ineff(DspSize9_TA_Ineff, :).SearchTime);

TP_prev_TA_Ineff_DspSize3_SEM = nanstd(TP_prev_TA_Ineff(DspSize3_TA_Ineff, :).SearchTime) / sqrt(length(DspSize3_TA_Ineff));
TP_prev_TA_Ineff_DspSize5_SEM = nanstd(TP_prev_TA_Ineff(DspSize5_TA_Ineff, :).SearchTime) / sqrt(length(DspSize5_TA_Ineff));
TP_prev_TA_Ineff_DspSize7_SEM = nanstd(TP_prev_TA_Ineff(DspSize7_TA_Ineff, :).SearchTime) / sqrt(length(DspSize7_TA_Ineff));
TP_prev_TA_Ineff_DspSize9_SEM = nanstd(TP_prev_TA_Ineff(DspSize9_TA_Ineff, :).SearchTime) / sqrt(length(DspSize9_TA_Ineff));


% Prepare data for plotting
display_sizes = [3, 5, 7, 9];

% Mean search times for Eff and Ineff for both TP_prev_TP and TP_prev_TA
mean_search_times_TP_Eff = [TP_prev_TP_Eff_DspSize3_SearchTime, TP_prev_TP_Eff_DspSize5_SearchTime, ...
                            TP_prev_TP_Eff_DspSize7_SearchTime, TP_prev_TP_Eff_DspSize9_SearchTime];
mean_search_times_TP_Ineff = [TP_prev_TP_Ineff_DspSize3_SearchTime, TP_prev_TP_Ineff_DspSize5_SearchTime, ...
                              TP_prev_TP_Ineff_DspSize7_SearchTime, TP_prev_TP_Ineff_DspSize9_SearchTime];

% SEM search times for Eff and Ineff for both TP_prev_TP and TP_prev_TA
sem_search_times_TP_Eff = [TP_prev_TP_Eff_DspSize3_SEM, TP_prev_TP_Eff_DspSize5_SEM, ...
                           TP_prev_TP_Eff_DspSize7_SEM, TP_prev_TP_Eff_DspSize9_SEM];
sem_search_times_TP_Ineff = [TP_prev_TP_Ineff_DspSize3_SEM, TP_prev_TP_Ineff_DspSize5_SEM, ...
                             TP_prev_TP_Ineff_DspSize7_SEM, TP_prev_TP_Ineff_DspSize9_SEM];


mean_search_times_TA_Eff = [TP_prev_TA_Eff_DspSize3_SearchTime, TP_prev_TA_Eff_DspSize5_SearchTime, ...
                            TP_prev_TA_Eff_DspSize7_SearchTime, TP_prev_TA_Eff_DspSize9_SearchTime];
mean_search_times_TA_Ineff = [TP_prev_TA_Ineff_DspSize3_SearchTime, TP_prev_TA_Ineff_DspSize5_SearchTime, ...
                              TP_prev_TA_Ineff_DspSize7_SearchTime, TP_prev_TA_Ineff_DspSize9_SearchTime];

% SEM search times for Eff and Ineff for both TP_prev_TP and TP_prev_TA
sem_search_times_TA_Eff = [TP_prev_TA_Eff_DspSize3_SEM, TP_prev_TA_Eff_DspSize5_SEM, ...
                           TP_prev_TA_Eff_DspSize7_SEM, TP_prev_TA_Eff_DspSize9_SEM];
sem_search_times_TA_Ineff = [TP_prev_TA_Ineff_DspSize3_SEM, TP_prev_TA_Ineff_DspSize5_SEM, ...
                             TP_prev_TA_Ineff_DspSize7_SEM, TP_prev_TA_Ineff_DspSize9_SEM];



% Plot with error bars for Eff and Ineff
figure;


% TP Prev TP Ineff (blue square dashed line)
errorbar(display_sizes, mean_search_times_TP_Ineff, sem_search_times_TP_Ineff, 'o-', 'LineWidth', 2, 'DisplayName', 'TP Prev TP Ineff');
hold on;



% TP Prev TA Ineff (red square dashed line)
errorbar(display_sizes, mean_search_times_TA_Ineff, sem_search_times_TA_Ineff, 's-', 'LineWidth', 2, 'DisplayName', 'TP Prev TA Ineff');
hold on;
set(gca,'FontSize',23); 


% Add labels and title
xlabel('Display Size');
ylabel('Mean Search Time');
title('Mean Search Time by Display Size for Inefficient TP prev TP and TP prev TA');
legend('show');
xlim([2,10]);

grid on;
hold off;

figure;
% TP Prev TP Eff (blue circle solid line)
errorbar(display_sizes, mean_search_times_TP_Eff, sem_search_times_TP_Eff, 'o-', 'LineWidth', 2, 'DisplayName', 'TP Prev TP Eff');
hold on;

% TP Prev TA Eff (red circle solid line)
errorbar(display_sizes, mean_search_times_TA_Eff, sem_search_times_TA_Eff, 's-', 'LineWidth', 2, 'DisplayName', 'TP Prev TA Eff');
hold on;

set(gca,'FontSize',23); 


% Add labels and title
xlabel('Display Size');
ylabel('Mean Search Time');
title('Mean Search Time by Display Size for Efficient TP prev TP and TP prev TA');
legend('show');
xlim([2,10]);

grid on;
hold off;

%% TA is Current trial - EFF/INEFF separated


% Separate Eff and Ineff conditions for TA_prev_TP
Eff_TP = find(TA_prev_TP.Search_Type == 1);
Ineff_TP = find(TA_prev_TP.Search_Type == 0);

TA_prev_TP_Eff = TA_prev_TP(Eff_TP, :);
TA_prev_TP_Ineff = TA_prev_TP(Ineff_TP, :);

% TA_TP display size separation for Eff
DspSize3_TP_Eff = find(TA_prev_TP_Eff.DispSize == 3);
DspSize5_TP_Eff = find(TA_prev_TP_Eff.DispSize == 5);
DspSize7_TP_Eff = find(TA_prev_TP_Eff.DispSize == 7);
DspSize9_TP_Eff = find(TA_prev_TP_Eff.DispSize == 9);

% Mean and SEM for Eff condition
TA_prev_TP_Eff_DspSize3_SearchTime = nanmean(TA_prev_TP_Eff(DspSize3_TP_Eff, :).SearchTime);
TA_prev_TP_Eff_DspSize5_SearchTime = nanmean(TA_prev_TP_Eff(DspSize5_TP_Eff, :).SearchTime);
TA_prev_TP_Eff_DspSize7_SearchTime = nanmean(TA_prev_TP_Eff(DspSize7_TP_Eff, :).SearchTime);
TA_prev_TP_Eff_DspSize9_SearchTime = nanmean(TA_prev_TP_Eff(DspSize9_TP_Eff, :).SearchTime);

TA_prev_TP_Eff_DspSize3_SEM = nanstd(TA_prev_TP_Eff(DspSize3_TP_Eff, :).SearchTime) / sqrt(length(DspSize3_TP_Eff));
TA_prev_TP_Eff_DspSize5_SEM = nanstd(TA_prev_TP_Eff(DspSize5_TP_Eff, :).SearchTime) / sqrt(length(DspSize5_TP_Eff));
TA_prev_TP_Eff_DspSize7_SEM = nanstd(TA_prev_TP_Eff(DspSize7_TP_Eff, :).SearchTime) / sqrt(length(DspSize7_TP_Eff));
TA_prev_TP_Eff_DspSize9_SEM = nanstd(TA_prev_TP_Eff(DspSize9_TP_Eff, :).SearchTime) / sqrt(length(DspSize9_TP_Eff));

% TA_TP display size separation for Ineff
DspSize3_TP_Ineff = find(TA_prev_TP_Ineff.DispSize == 3);
DspSize5_TP_Ineff = find(TA_prev_TP_Ineff.DispSize == 5);
DspSize7_TP_Ineff = find(TA_prev_TP_Ineff.DispSize == 7);
DspSize9_TP_Ineff = find(TA_prev_TP_Ineff.DispSize == 9);

% Mean and SEM for Ineff condition
TA_prev_TP_Ineff_DspSize3_SearchTime = nanmean(TA_prev_TP_Ineff(DspSize3_TP_Ineff, :).SearchTime);
TA_prev_TP_Ineff_DspSize5_SearchTime = nanmean(TA_prev_TP_Ineff(DspSize5_TP_Ineff, :).SearchTime);
TA_prev_TP_Ineff_DspSize7_SearchTime = nanmean(TA_prev_TP_Ineff(DspSize7_TP_Ineff, :).SearchTime);
TA_prev_TP_Ineff_DspSize9_SearchTime = nanmean(TA_prev_TP_Ineff(DspSize9_TP_Ineff, :).SearchTime);

TA_prev_TP_Ineff_DspSize3_SEM = nanstd(TA_prev_TP_Ineff(DspSize3_TP_Ineff, :).SearchTime) / sqrt(length(DspSize3_TP_Ineff));
TA_prev_TP_Ineff_DspSize5_SEM = nanstd(TA_prev_TP_Ineff(DspSize5_TP_Ineff, :).SearchTime) / sqrt(length(DspSize5_TP_Ineff));
TA_prev_TP_Ineff_DspSize7_SEM = nanstd(TA_prev_TP_Ineff(DspSize7_TP_Ineff, :).SearchTime) / sqrt(length(DspSize7_TP_Ineff));
TA_prev_TP_Ineff_DspSize9_SEM = nanstd(TA_prev_TP_Ineff(DspSize9_TP_Ineff, :).SearchTime) / sqrt(length(DspSize9_TP_Ineff));

% Separate Eff and Ineff conditions for TA_prev_TA
Eff_TA = find(TA_prev_TA.Search_Type == 1);
Ineff_TA = find(TA_prev_TA.Search_Type == 0);

TA_prev_TA_Eff = TA_prev_TA(Eff_TA, :);
TA_prev_TA_Ineff = TA_prev_TA(Ineff_TA, :);

% TA_TA display size separation for Eff
DspSize3_TA_Eff = find(TA_prev_TA_Eff.DispSize == 3);
DspSize5_TA_Eff = find(TA_prev_TA_Eff.DispSize == 5);
DspSize7_TA_Eff = find(TA_prev_TA_Eff.DispSize == 7);
DspSize9_TA_Eff = find(TA_prev_TA_Eff.DispSize == 9);

% Mean and SEM for Eff condition
TA_prev_TA_Eff_DspSize3_SearchTime = nanmean(TA_prev_TA_Eff(DspSize3_TA_Eff, :).SearchTime);
TA_prev_TA_Eff_DspSize5_SearchTime = nanmean(TA_prev_TA_Eff(DspSize5_TA_Eff, :).SearchTime);
TA_prev_TA_Eff_DspSize7_SearchTime = nanmean(TA_prev_TA_Eff(DspSize7_TA_Eff, :).SearchTime);
TA_prev_TA_Eff_DspSize9_SearchTime = nanmean(TA_prev_TA_Eff(DspSize9_TA_Eff, :).SearchTime);

TA_prev_TA_Eff_DspSize3_SEM = nanstd(TA_prev_TA_Eff(DspSize3_TA_Eff, :).SearchTime) / sqrt(length(DspSize3_TA_Eff));
TA_prev_TA_Eff_DspSize5_SEM = nanstd(TA_prev_TA_Eff(DspSize5_TA_Eff, :).SearchTime) / sqrt(length(DspSize5_TA_Eff));
TA_prev_TA_Eff_DspSize7_SEM = nanstd(TA_prev_TA_Eff(DspSize7_TA_Eff, :).SearchTime) / sqrt(length(DspSize7_TA_Eff));
TA_prev_TA_Eff_DspSize9_SEM = nanstd(TA_prev_TA_Eff(DspSize9_TA_Eff, :).SearchTime) / sqrt(length(DspSize9_TA_Eff));

% TA_TA display size separation for Ineff
DspSize3_TA_Ineff = find(TA_prev_TA_Ineff.DispSize == 3);
DspSize5_TA_Ineff = find(TA_prev_TA_Ineff.DispSize == 5);
DspSize7_TA_Ineff = find(TA_prev_TA_Ineff.DispSize == 7);
DspSize9_TA_Ineff = find(TA_prev_TA_Ineff.DispSize == 9);

% Mean and SEM for Ineff condition
TA_prev_TA_Ineff_DspSize3_SearchTime = nanmean(TA_prev_TA_Ineff(DspSize3_TA_Ineff, :).SearchTime);
TA_prev_TA_Ineff_DspSize5_SearchTime = nanmean(TA_prev_TA_Ineff(DspSize5_TA_Ineff, :).SearchTime);
TA_prev_TA_Ineff_DspSize7_SearchTime = nanmean(TA_prev_TA_Ineff(DspSize7_TA_Ineff, :).SearchTime);
TA_prev_TA_Ineff_DspSize9_SearchTime = nanmean(TA_prev_TA_Ineff(DspSize9_TA_Ineff, :).SearchTime);

TA_prev_TA_Ineff_DspSize3_SEM = nanstd(TA_prev_TA_Ineff(DspSize3_TA_Ineff, :).SearchTime) / sqrt(length(DspSize3_TA_Ineff));
TA_prev_TA_Ineff_DspSize5_SEM = nanstd(TA_prev_TA_Ineff(DspSize5_TA_Ineff, :).SearchTime) / sqrt(length(DspSize5_TA_Ineff));
TA_prev_TA_Ineff_DspSize7_SEM = nanstd(TA_prev_TA_Ineff(DspSize7_TA_Ineff, :).SearchTime) / sqrt(length(DspSize7_TA_Ineff));
TA_prev_TA_Ineff_DspSize9_SEM = nanstd(TA_prev_TA_Ineff(DspSize9_TA_Ineff, :).SearchTime) / sqrt(length(DspSize9_TA_Ineff));


% Prepare data for plotting
display_sizes = [3, 5, 7, 9];

% Mean search times for Eff and Ineff for both TA_prev_TA and TA_prev_TP
mean_search_times_TP_Eff = [TA_prev_TP_Eff_DspSize3_SearchTime, TA_prev_TP_Eff_DspSize5_SearchTime, ...
                            TA_prev_TP_Eff_DspSize7_SearchTime, TA_prev_TP_Eff_DspSize9_SearchTime];
mean_search_times_TP_Ineff = [TA_prev_TP_Ineff_DspSize3_SearchTime, TA_prev_TP_Ineff_DspSize5_SearchTime, ...
                              TA_prev_TP_Ineff_DspSize7_SearchTime, TA_prev_TP_Ineff_DspSize9_SearchTime];

mean_search_times_TA_Eff = [TA_prev_TA_Eff_DspSize3_SearchTime, TA_prev_TA_Eff_DspSize5_SearchTime, ...
                            TA_prev_TA_Eff_DspSize7_SearchTime, TA_prev_TA_Eff_DspSize9_SearchTime];
mean_search_times_TA_Ineff = [TA_prev_TA_Ineff_DspSize3_SearchTime, TA_prev_TA_Ineff_DspSize5_SearchTime, ...
                              TA_prev_TA_Ineff_DspSize7_SearchTime, TA_prev_TA_Ineff_DspSize9_SearchTime];

% SEM search times for Eff and Ineff for both TA_prev_TA and TA_prev_TP
sem_search_times_TP_Eff = [TA_prev_TP_Eff_DspSize3_SEM, TA_prev_TP_Eff_DspSize5_SEM, ...
                           TA_prev_TP_Eff_DspSize7_SEM, TA_prev_TP_Eff_DspSize9_SEM];
sem_search_times_TP_Ineff = [TA_prev_TP_Ineff_DspSize3_SEM, TA_prev_TP_Ineff_DspSize5_SEM, ...
                             TA_prev_TP_Ineff_DspSize7_SEM, TA_prev_TP_Ineff_DspSize9_SEM];

sem_search_times_TA_Eff = [TA_prev_TA_Eff_DspSize3_SEM, TA_prev_TA_Eff_DspSize5_SEM, ...
                           TA_prev_TA_Eff_DspSize7_SEM, TA_prev_TA_Eff_DspSize9_SEM];
sem_search_times_TA_Ineff = [TA_prev_TA_Ineff_DspSize3_SEM, TA_prev_TA_Ineff_DspSize5_SEM, ...
                             TA_prev_TA_Ineff_DspSize7_SEM, TA_prev_TA_Ineff_DspSize9_SEM];


% Plot with error bars for Eff and Ineff
figure;


% TA Prev TP Ineff (blue square dashed line)
errorbar(display_sizes, mean_search_times_TP_Ineff, sem_search_times_TP_Ineff, 'o-', 'LineWidth', 2, 'DisplayName', 'TA Prev TP Ineff');
hold on;


% TA Prev TA Ineff (red square dashed line)
errorbar(display_sizes, mean_search_times_TA_Ineff, sem_search_times_TA_Ineff, 's-', 'LineWidth', 2, 'DisplayName', 'TA Prev TA Ineff');
hold on;

set(gca,'FontSize',23); 

% Add labels and title
xlabel('Display Size');
ylabel('Mean Search Time');
title('Mean Search Time by Display Size for Inefficient TA prev TP and TA prev TA');
legend('show');
xlim([2,10]);

grid on;
hold off;


figure;
% TA Prev TP Eff (red circle solid line)
errorbar(display_sizes, mean_search_times_TP_Eff, sem_search_times_TP_Eff, 'o-', 'LineWidth', 2, 'DisplayName', 'TA Prev TP Eff');

hold on;

% TA Prev TA Eff (blue circle solid line)
errorbar(display_sizes, mean_search_times_TA_Eff, sem_search_times_TA_Eff, 's-', 'LineWidth', 2, 'DisplayName', 'TA Prev TA Eff');

hold on;

set(gca,'FontSize',23); 

% Add labels and title
xlabel('Display Size');
ylabel('Mean Search Time');
title('Mean Search Time by Display Size for Efficient TA prev TP and TA prev TA');
legend('show');
xlim([2,10]);

grid on;
hold off;


%% %% Current trial TP - EFF/INEFF separated by Subject (Monkeys 3 and 4)

% Separate Eff and Ineff conditions for each monkey
% Monkey 3
Eff_TP_Monkey3 = find(TP_prev_TP_Final.Subject == 3 & TP_prev_TP_Final.Search_Type == 1);
Ineff_TP_Monkey3 = find(TP_prev_TP_Final.Subject == 3 & TP_prev_TP_Final.Search_Type == 0);

TP_prev_TP_Final_Eff_Monkey3 = TP_prev_TP_Final(Eff_TP_Monkey3, :);
TP_prev_TP_Final_Ineff_Monkey3 = TP_prev_TP_Final(Ineff_TP_Monkey3, :);

% Monkey 4
Eff_TP_Monkey4 = find(TP_prev_TP_Final.Subject == 4 & TP_prev_TP_Final.Search_Type == 1);
Ineff_TP_Monkey4 = find(TP_prev_TP_Final.Subject == 4 & TP_prev_TP_Final.Search_Type == 0);

TP_prev_TP_Final_Eff_Monkey4 = TP_prev_TP_Final(Eff_TP_Monkey4, :);
TP_prev_TP_Final_Ineff_Monkey4 = TP_prev_TP_Final(Ineff_TP_Monkey4, :);

% Separate Eff and Ineff conditions for TP_prev_TA_Final for each monkey
% Monkey 3
Eff_TA_Monkey3 = find(TP_prev_TA_Final.Subject == 3 & TP_prev_TA_Final.Search_Type == 1);
Ineff_TA_Monkey3 = find(TP_prev_TA_Final.Subject == 3 & TP_prev_TA_Final.Search_Type == 0);

TP_prev_TA_Final_Eff_Monkey3 = TP_prev_TA_Final(Eff_TA_Monkey3, :);
TP_prev_TA_Final_Ineff_Monkey3 = TP_prev_TA_Final(Ineff_TA_Monkey3, :);

% Monkey 4
Eff_TA_Monkey4 = find(TP_prev_TA_Final.Subject == 4 & TP_prev_TA_Final.Search_Type == 1);
Ineff_TA_Monkey4 = find(TP_prev_TA_Final.Subject == 4 & TP_prev_TA_Final.Search_Type == 0);

TP_prev_TA_Final_Eff_Monkey4 = TP_prev_TA_Final(Eff_TA_Monkey4, :);
TP_prev_TA_Final_Ineff_Monkey4 = TP_prev_TA_Final(Ineff_TA_Monkey4, :);

% Display size separation for Eff and Ineff for Monkey 3 and Monkey 4

% TP_TP display size separation for Eff for Monkey 3
DspSize3_TP_Eff_Monkey3 = find(TP_prev_TP_Final_Eff_Monkey3.DispSize == 3);
DspSize5_TP_Eff_Monkey3 = find(TP_prev_TP_Final_Eff_Monkey3.DispSize == 5);
DspSize7_TP_Eff_Monkey3 = find(TP_prev_TP_Final_Eff_Monkey3.DispSize == 7);
DspSize9_TP_Eff_Monkey3 = find(TP_prev_TP_Final_Eff_Monkey3.DispSize == 9);

Monkey3_TP_Eff_SearchTime_3 = nanmean(TP_prev_TP_Final_Eff_Monkey3(DspSize3_TP_Eff_Monkey3, :).SearchTime);
Monkey3_TP_Eff_SearchTime_5 = nanmean(TP_prev_TP_Final_Eff_Monkey3(DspSize5_TP_Eff_Monkey3, :).SearchTime);
Monkey3_TP_Eff_SearchTime_7 = nanmean(TP_prev_TP_Final_Eff_Monkey3(DspSize7_TP_Eff_Monkey3, :).SearchTime);
Monkey3_TP_Eff_SearchTime_9 = nanmean(TP_prev_TP_Final_Eff_Monkey3(DspSize9_TP_Eff_Monkey3, :).SearchTime);

Monkey3_TP_Eff_SEM_3 = nanstd(TP_prev_TP_Final_Eff_Monkey3(DspSize3_TP_Eff_Monkey3, :).SearchTime) / sqrt(length(DspSize3_TP_Eff_Monkey3));
Monkey3_TP_Eff_SEM_5 = nanstd(TP_prev_TP_Final_Eff_Monkey3(DspSize5_TP_Eff_Monkey3, :).SearchTime) / sqrt(length(DspSize5_TP_Eff_Monkey3));
Monkey3_TP_Eff_SEM_7 = nanstd(TP_prev_TP_Final_Eff_Monkey3(DspSize7_TP_Eff_Monkey3, :).SearchTime) / sqrt(length(DspSize7_TP_Eff_Monkey3));
Monkey3_TP_Eff_SEM_9 = nanstd(TP_prev_TP_Final_Eff_Monkey3(DspSize9_TP_Eff_Monkey3, :).SearchTime) / sqrt(length(DspSize9_TP_Eff_Monkey3));

% Repeat the same for Ineff for Monkey 3
DspSize3_TP_Ineff_Monkey3 = find(TP_prev_TP_Final_Ineff_Monkey3.DispSize == 3);
DspSize5_TP_Ineff_Monkey3 = find(TP_prev_TP_Final_Ineff_Monkey3.DispSize == 5);
DspSize7_TP_Ineff_Monkey3 = find(TP_prev_TP_Final_Ineff_Monkey3.DispSize == 7);
DspSize9_TP_Ineff_Monkey3 = find(TP_prev_TP_Final_Ineff_Monkey3.DispSize == 9);

Monkey3_TP_Ineff_SearchTime_3 = nanmean(TP_prev_TP_Final_Ineff_Monkey3(DspSize3_TP_Ineff_Monkey3, :).SearchTime);
Monkey3_TP_Ineff_SearchTime_5 = nanmean(TP_prev_TP_Final_Ineff_Monkey3(DspSize5_TP_Ineff_Monkey3, :).SearchTime);
Monkey3_TP_Ineff_SearchTime_7 = nanmean(TP_prev_TP_Final_Ineff_Monkey3(DspSize7_TP_Ineff_Monkey3, :).SearchTime);
Monkey3_TP_Ineff_SearchTime_9 = nanmean(TP_prev_TP_Final_Ineff_Monkey3(DspSize9_TP_Ineff_Monkey3, :).SearchTime);

Monkey3_TP_Ineff_SEM_3 = nanstd(TP_prev_TP_Final_Ineff_Monkey3(DspSize3_TP_Ineff_Monkey3, :).SearchTime) / sqrt(length(DspSize3_TP_Ineff_Monkey3));
Monkey3_TP_Ineff_SEM_5 = nanstd(TP_prev_TP_Final_Ineff_Monkey3(DspSize5_TP_Ineff_Monkey3, :).SearchTime) / sqrt(length(DspSize5_TP_Ineff_Monkey3));
Monkey3_TP_Ineff_SEM_7 = nanstd(TP_prev_TP_Final_Ineff_Monkey3(DspSize7_TP_Ineff_Monkey3, :).SearchTime) / sqrt(length(DspSize7_TP_Ineff_Monkey3));
Monkey3_TP_Ineff_SEM_9 = nanstd(TP_prev_TP_Final_Ineff_Monkey3(DspSize9_TP_Ineff_Monkey3, :).SearchTime) / sqrt(length(DspSize9_TP_Ineff_Monkey3));

% Display size separation for Eff and Ineff for Monkey4

% TP_TP display size separation for Eff for Monkey 4
DspSize3_TP_Eff_Monkey4 = find(TP_prev_TP_Final_Eff_Monkey4.DispSize == 3);
DspSize5_TP_Eff_Monkey4 = find(TP_prev_TP_Final_Eff_Monkey4.DispSize == 5);
DspSize7_TP_Eff_Monkey4 = find(TP_prev_TP_Final_Eff_Monkey4.DispSize == 7);
DspSize9_TP_Eff_Monkey4 = find(TP_prev_TP_Final_Eff_Monkey4.DispSize == 9);

Monkey4_TP_Eff_SearchTime_3 = nanmean(TP_prev_TP_Final_Eff_Monkey4(DspSize3_TP_Eff_Monkey4, :).SearchTime);
Monkey4_TP_Eff_SearchTime_5 = nanmean(TP_prev_TP_Final_Eff_Monkey4(DspSize5_TP_Eff_Monkey4, :).SearchTime);
Monkey4_TP_Eff_SearchTime_7 = nanmean(TP_prev_TP_Final_Eff_Monkey4(DspSize7_TP_Eff_Monkey4, :).SearchTime);
Monkey4_TP_Eff_SearchTime_9 = nanmean(TP_prev_TP_Final_Eff_Monkey4(DspSize9_TP_Eff_Monkey4, :).SearchTime);

Monkey4_TP_Eff_SEM_3 = nanstd(TP_prev_TP_Final_Eff_Monkey4(DspSize3_TP_Eff_Monkey4, :).SearchTime) / sqrt(length(DspSize3_TP_Eff_Monkey4));
Monkey4_TP_Eff_SEM_5 = nanstd(TP_prev_TP_Final_Eff_Monkey4(DspSize5_TP_Eff_Monkey4, :).SearchTime) / sqrt(length(DspSize5_TP_Eff_Monkey4));
Monkey4_TP_Eff_SEM_7 = nanstd(TP_prev_TP_Final_Eff_Monkey4(DspSize7_TP_Eff_Monkey4, :).SearchTime) / sqrt(length(DspSize7_TP_Eff_Monkey4));
Monkey4_TP_Eff_SEM_9 = nanstd(TP_prev_TP_Final_Eff_Monkey4(DspSize9_TP_Eff_Monkey4, :).SearchTime) / sqrt(length(DspSize9_TP_Eff_Monkey4));

% Repeat the same for Ineff for Monkey 4
DspSize3_TP_Ineff_Monkey4 = find(TP_prev_TP_Final_Ineff_Monkey4.DispSize == 3);
DspSize5_TP_Ineff_Monkey4 = find(TP_prev_TP_Final_Ineff_Monkey4.DispSize == 5);
DspSize7_TP_Ineff_Monkey4 = find(TP_prev_TP_Final_Ineff_Monkey4.DispSize == 7);
DspSize9_TP_Ineff_Monkey4 = find(TP_prev_TP_Final_Ineff_Monkey4.DispSize == 9);

Monkey4_TP_Ineff_SearchTime_3 = nanmean(TP_prev_TP_Final_Ineff_Monkey4(DspSize3_TP_Ineff_Monkey4, :).SearchTime);
Monkey4_TP_Ineff_SearchTime_5 = nanmean(TP_prev_TP_Final_Ineff_Monkey4(DspSize5_TP_Ineff_Monkey4, :).SearchTime);
Monkey4_TP_Ineff_SearchTime_7 = nanmean(TP_prev_TP_Final_Ineff_Monkey4(DspSize7_TP_Ineff_Monkey4, :).SearchTime);
Monkey4_TP_Ineff_SearchTime_9 = nanmean(TP_prev_TP_Final_Ineff_Monkey4(DspSize9_TP_Ineff_Monkey4, :).SearchTime);

Monkey4_TP_Ineff_SEM_3 = nanstd(TP_prev_TP_Final_Ineff_Monkey4(DspSize3_TP_Ineff_Monkey4, :).SearchTime) / sqrt(length(DspSize3_TP_Ineff_Monkey4));
Monkey4_TP_Ineff_SEM_5 = nanstd(TP_prev_TP_Final_Ineff_Monkey4(DspSize5_TP_Ineff_Monkey4, :).SearchTime) / sqrt(length(DspSize5_TP_Ineff_Monkey4));
Monkey4_TP_Ineff_SEM_7 = nanstd(TP_prev_TP_Final_Ineff_Monkey4(DspSize7_TP_Ineff_Monkey4, :).SearchTime) / sqrt(length(DspSize7_TP_Ineff_Monkey4));
Monkey4_TP_Ineff_SEM_9 = nanstd(TP_prev_TP_Final_Ineff_Monkey4(DspSize9_TP_Ineff_Monkey4, :).SearchTime) / sqrt(length(DspSize9_TP_Ineff_Monkey4));

% Display size separation for Eff and Ineff for Monkey 3 and Monkey 4 for TP Prev TA

% TP_TA display size separation for Eff for Monkey 3
DspSize3_TA_Eff_Monkey3 = find(TP_prev_TA_Final_Eff_Monkey3.DispSize == 3);
DspSize5_TA_Eff_Monkey3 = find(TP_prev_TA_Final_Eff_Monkey3.DispSize == 5);
DspSize7_TA_Eff_Monkey3 = find(TP_prev_TA_Final_Eff_Monkey3.DispSize == 7);
DspSize9_TA_Eff_Monkey3 = find(TP_prev_TA_Final_Eff_Monkey3.DispSize == 9);

Monkey3_TA_Eff_SearchTime_3 = nanmean(TP_prev_TA_Final_Eff_Monkey3(DspSize3_TA_Eff_Monkey3, :).SearchTime);
Monkey3_TA_Eff_SearchTime_5 = nanmean(TP_prev_TA_Final_Eff_Monkey3(DspSize5_TA_Eff_Monkey3, :).SearchTime);
Monkey3_TA_Eff_SearchTime_7 = nanmean(TP_prev_TA_Final_Eff_Monkey3(DspSize7_TA_Eff_Monkey3, :).SearchTime);
Monkey3_TA_Eff_SearchTime_9 = nanmean(TP_prev_TA_Final_Eff_Monkey3(DspSize9_TA_Eff_Monkey3, :).SearchTime);

Monkey3_TA_Eff_SEM_3 = nanstd(TP_prev_TA_Final_Eff_Monkey3(DspSize3_TA_Eff_Monkey3, :).SearchTime) / sqrt(length(DspSize3_TA_Eff_Monkey3));
Monkey3_TA_Eff_SEM_5 = nanstd(TP_prev_TA_Final_Eff_Monkey3(DspSize5_TA_Eff_Monkey3, :).SearchTime) / sqrt(length(DspSize5_TA_Eff_Monkey3));
Monkey3_TA_Eff_SEM_7 = nanstd(TP_prev_TA_Final_Eff_Monkey3(DspSize7_TA_Eff_Monkey3, :).SearchTime) / sqrt(length(DspSize7_TA_Eff_Monkey3));
Monkey3_TA_Eff_SEM_9 = nanstd(TP_prev_TA_Final_Eff_Monkey3(DspSize9_TA_Eff_Monkey3, :).SearchTime) / sqrt(length(DspSize9_TA_Eff_Monkey3));

% Repeat the same for Ineff for Monkey 3
DspSize3_TA_Ineff_Monkey3 = find(TP_prev_TA_Final_Ineff_Monkey3.DispSize == 3);
DspSize5_TA_Ineff_Monkey3 = find(TP_prev_TA_Final_Ineff_Monkey3.DispSize == 5);
DspSize7_TA_Ineff_Monkey3 = find(TP_prev_TA_Final_Ineff_Monkey3.DispSize == 7);
DspSize9_TA_Ineff_Monkey3 = find(TP_prev_TA_Final_Ineff_Monkey3.DispSize == 9);

Monkey3_TA_Ineff_SearchTime_3 = nanmean(TP_prev_TA_Final_Ineff_Monkey3(DspSize3_TA_Ineff_Monkey3, :).SearchTime);
Monkey3_TA_Ineff_SearchTime_5 = nanmean(TP_prev_TA_Final_Ineff_Monkey3(DspSize5_TA_Ineff_Monkey3, :).SearchTime);
Monkey3_TA_Ineff_SearchTime_7 = nanmean(TP_prev_TA_Final_Ineff_Monkey3(DspSize7_TA_Ineff_Monkey3, :).SearchTime);
Monkey3_TA_Ineff_SearchTime_9 = nanmean(TP_prev_TA_Final_Ineff_Monkey3(DspSize9_TA_Ineff_Monkey3, :).SearchTime);

Monkey3_TA_Ineff_SEM_3 = nanstd(TP_prev_TA_Final_Ineff_Monkey3(DspSize3_TA_Ineff_Monkey3, :).SearchTime) / sqrt(length(DspSize3_TA_Ineff_Monkey3));
Monkey3_TA_Ineff_SEM_5 = nanstd(TP_prev_TA_Final_Ineff_Monkey3(DspSize5_TA_Ineff_Monkey3, :).SearchTime) / sqrt(length(DspSize5_TA_Ineff_Monkey3));
Monkey3_TA_Ineff_SEM_7 = nanstd(TP_prev_TA_Final_Ineff_Monkey3(DspSize7_TA_Ineff_Monkey3, :).SearchTime) / sqrt(length(DspSize7_TA_Ineff_Monkey3));
Monkey3_TA_Ineff_SEM_9 = nanstd(TP_prev_TA_Final_Ineff_Monkey3(DspSize9_TA_Ineff_Monkey3, :).SearchTime) / sqrt(length(DspSize9_TA_Ineff_Monkey3));

% Display size separation for Eff and Ineff for Monkey 4

% TP_TA display size separation for Eff for Monkey 4
DspSize3_TA_Eff_Monkey4 = find(TP_prev_TA_Final_Eff_Monkey4.DispSize == 3);
DspSize5_TA_Eff_Monkey4 = find(TP_prev_TA_Final_Eff_Monkey4.DispSize == 5);
DspSize7_TA_Eff_Monkey4 = find(TP_prev_TA_Final_Eff_Monkey4.DispSize == 7);
DspSize9_TA_Eff_Monkey4 = find(TP_prev_TA_Final_Eff_Monkey4.DispSize == 9);

Monkey4_TA_Eff_SearchTime_3 = nanmean(TP_prev_TA_Final_Eff_Monkey4(DspSize3_TA_Eff_Monkey4, :).SearchTime);
Monkey4_TA_Eff_SearchTime_5 = nanmean(TP_prev_TA_Final_Eff_Monkey4(DspSize5_TA_Eff_Monkey4, :).SearchTime);
Monkey4_TA_Eff_SearchTime_7 = nanmean(TP_prev_TA_Final_Eff_Monkey4(DspSize7_TA_Eff_Monkey4, :).SearchTime);
Monkey4_TA_Eff_SearchTime_9 = nanmean(TP_prev_TA_Final_Eff_Monkey4(DspSize9_TA_Eff_Monkey4, :).SearchTime);

Monkey4_TA_Eff_SEM_3 = nanstd(TP_prev_TA_Final_Eff_Monkey4(DspSize3_TA_Eff_Monkey4, :).SearchTime) / sqrt(length(DspSize3_TA_Eff_Monkey4));
Monkey4_TA_Eff_SEM_5 = nanstd(TP_prev_TA_Final_Eff_Monkey4(DspSize5_TA_Eff_Monkey4, :).SearchTime) / sqrt(length(DspSize5_TA_Eff_Monkey4));
Monkey4_TA_Eff_SEM_7 = nanstd(TP_prev_TA_Final_Eff_Monkey4(DspSize7_TA_Eff_Monkey4, :).SearchTime) / sqrt(length(DspSize7_TA_Eff_Monkey4));
Monkey4_TA_Eff_SEM_9 = nanstd(TP_prev_TA_Final_Eff_Monkey4(DspSize9_TA_Eff_Monkey4, :).SearchTime) / sqrt(length(DspSize9_TA_Eff_Monkey4));

% Repeat the same for Ineff for Monkey 4
DspSize3_TA_Ineff_Monkey4 = find(TP_prev_TA_Final_Ineff_Monkey4.DispSize == 3);
DspSize5_TA_Ineff_Monkey4 = find(TP_prev_TA_Final_Ineff_Monkey4.DispSize == 5);
DspSize7_TA_Ineff_Monkey4 = find(TP_prev_TA_Final_Ineff_Monkey4.DispSize == 7);
DspSize9_TA_Ineff_Monkey4 = find(TP_prev_TA_Final_Ineff_Monkey4.DispSize == 9);

Monkey4_TA_Ineff_SearchTime_3 = nanmean(TP_prev_TA_Final_Ineff_Monkey4(DspSize3_TA_Ineff_Monkey4, :).SearchTime);
Monkey4_TA_Ineff_SearchTime_5 = nanmean(TP_prev_TA_Final_Ineff_Monkey4(DspSize5_TA_Ineff_Monkey4, :).SearchTime);
Monkey4_TA_Ineff_SearchTime_7 = nanmean(TP_prev_TA_Final_Ineff_Monkey4(DspSize7_TA_Ineff_Monkey4, :).SearchTime);
Monkey4_TA_Ineff_SearchTime_9 = nanmean(TP_prev_TA_Final_Ineff_Monkey4(DspSize9_TA_Ineff_Monkey4, :).SearchTime);

Monkey4_TA_Ineff_SEM_3 = nanstd(TP_prev_TA_Final_Ineff_Monkey4(DspSize3_TA_Ineff_Monkey4, :).SearchTime) / sqrt(length(DspSize3_TA_Ineff_Monkey4));
Monkey4_TA_Ineff_SEM_5 = nanstd(TP_prev_TA_Final_Ineff_Monkey4(DspSize5_TA_Ineff_Monkey4, :).SearchTime) / sqrt(length(DspSize5_TA_Ineff_Monkey4));
Monkey4_TA_Ineff_SEM_7 = nanstd(TP_prev_TA_Final_Ineff_Monkey4(DspSize7_TA_Ineff_Monkey4, :).SearchTime) / sqrt(length(DspSize7_TA_Ineff_Monkey4));
Monkey4_TA_Ineff_SEM_9 = nanstd(TP_prev_TA_Final_Ineff_Monkey4(DspSize9_TA_Ineff_Monkey4, :).SearchTime) / sqrt(length(DspSize9_TA_Ineff_Monkey4));

% Prepare data for plotting
display_sizes = [3, 5, 7, 9];

% Plot with error bars for Ineff
figure;

% TP Prev TP Ineff 
errorbar(display_sizes, [Monkey3_TP_Ineff_SearchTime_3, Monkey3_TP_Ineff_SearchTime_5, Monkey3_TP_Ineff_SearchTime_7, Monkey3_TP_Ineff_SearchTime_9], ...
    [Monkey3_TP_Ineff_SEM_3, Monkey3_TP_Ineff_SEM_5, Monkey3_TP_Ineff_SEM_7, Monkey3_TP_Ineff_SEM_9], 'o-', 'LineWidth', 2, 'DisplayName', 'TP Prev TP Ineff', 'Color', 'r');
hold on;
errorbar(display_sizes, [Monkey3_TA_Ineff_SearchTime_3, Monkey3_TA_Ineff_SearchTime_5, Monkey3_TA_Ineff_SearchTime_7, Monkey3_TA_Ineff_SearchTime_9], ...
    [Monkey3_TA_Ineff_SEM_3, Monkey3_TA_Ineff_SEM_5, Monkey3_TA_Ineff_SEM_7, Monkey3_TA_Ineff_SEM_9], 'o-', 'LineWidth', 2, 'DisplayName', 'TP Prev TA Ineff', 'Color', 'b');
xlabel('Display Size');
ylabel('Mean Search Time');
title('Mean Search Time by Display Size for Ineff (Monkeys 3)');
legend('show');
grid on;
hold off;

% Repeat similar plotting for Eff conditions
figure;

% TP Prev TP Eff (Monkey 3 in red, Monkey 4 in blue)
errorbar(display_sizes, [Monkey3_TP_Eff_SearchTime_3, Monkey3_TP_Eff_SearchTime_5, Monkey3_TP_Eff_SearchTime_7, Monkey3_TP_Eff_SearchTime_9], ...
    [Monkey3_TP_Eff_SEM_3, Monkey3_TP_Eff_SEM_5, Monkey3_TP_Eff_SEM_7, Monkey3_TP_Eff_SEM_9], 'o-', 'LineWidth', 2, 'DisplayName', 'TP Prev TP Eff', 'Color', 'r');
hold on;
errorbar(display_sizes, [Monkey3_TA_Eff_SearchTime_3, Monkey3_TA_Eff_SearchTime_5, Monkey3_TA_Eff_SearchTime_7, Monkey3_TA_Eff_SearchTime_9], ...
    [Monkey3_TA_Eff_SEM_3, Monkey3_TA_Eff_SEM_5, Monkey3_TA_Eff_SEM_7, Monkey3_TA_Eff_SEM_9], 'o-', 'LineWidth', 2, 'DisplayName', 'TP Prev TA Eff', 'Color', 'b');
xlabel('Display Size');
ylabel('Mean Search Time');
title('Mean Search Time by Display Size for Eff (Monkeys 3)');
legend('show');
grid on;
hold off;

% Repeat similar plotting for TP Prev TA Eff and Ineff conditions
figure;

% TP Prev TA Ineff (Monkey 3 in red, Monkey 4 in blue)
errorbar(display_sizes, [Monkey4_TP_Ineff_SearchTime_3, Monkey4_TP_Ineff_SearchTime_5, Monkey4_TP_Ineff_SearchTime_7, Monkey4_TP_Ineff_SearchTime_9], ...
    [Monkey4_TP_Ineff_SEM_3, Monkey4_TP_Ineff_SEM_5, Monkey4_TP_Ineff_SEM_7, Monkey4_TP_Ineff_SEM_9], 's-', 'LineWidth', 2, 'DisplayName', 'TP Prev TP Ineff', 'Color', 'r');
hold on;
errorbar(display_sizes, [Monkey4_TA_Ineff_SearchTime_3, Monkey4_TA_Ineff_SearchTime_5, Monkey4_TA_Ineff_SearchTime_7, Monkey4_TA_Ineff_SearchTime_9], ...
    [Monkey4_TA_Ineff_SEM_3, Monkey4_TA_Ineff_SEM_5, Monkey4_TA_Ineff_SEM_7, Monkey4_TA_Ineff_SEM_9], 's-', 'LineWidth', 2, 'DisplayName', 'TP Prev TA Ineff', 'Color', 'b');
xlabel('Display Size');
ylabel('Mean Search Time');
title('Mean Search Time by Display Size for Ineff (Monkey 4)');
legend('show');
grid on;
hold off;

figure;

% TP Prev TA Eff (Monkey 3 in red, Monkey 4 in blue)
errorbar(display_sizes, [Monkey4_TP_Eff_SearchTime_3, Monkey4_TP_Eff_SearchTime_5, Monkey4_TP_Eff_SearchTime_7, Monkey4_TP_Eff_SearchTime_9], ...
    [Monkey4_TP_Eff_SEM_3, Monkey4_TP_Eff_SEM_5, Monkey4_TP_Eff_SEM_7, Monkey4_TP_Eff_SEM_9], 's-', 'LineWidth', 2, 'DisplayName', 'TP Prev TP Eff', 'Color', 'r');
hold on;
errorbar(display_sizes, [Monkey4_TA_Eff_SearchTime_3, Monkey4_TA_Eff_SearchTime_5, Monkey4_TA_Eff_SearchTime_7, Monkey4_TA_Eff_SearchTime_9], ...
    [Monkey4_TA_Eff_SEM_3, Monkey4_TA_Eff_SEM_5, Monkey4_TA_Eff_SEM_7, Monkey4_TA_Eff_SEM_9], 's-', 'LineWidth', 2, 'DisplayName', 'TP Prev TA Eff', 'Color', 'b');
xlabel('Display Size');
ylabel('Mean Search Time');
title('Mean Search Time by Display Size for Eff (Monkey 4)');
legend('show');
grid on;
hold off;




%% Current trial TA - EFF/INEFF separated by Subject (Monkeys 3 and 4)

% Separate Eff and Ineff conditions for each monkey for TA_prev_TP_Final
% Monkey 3
Eff_TP_Monkey3 = find(TA_prev_TP_Final.Subject == 3 & TA_prev_TP_Final.Search_Type == 1);
Ineff_TP_Monkey3 = find(TA_prev_TP_Final.Subject == 3 & TA_prev_TP_Final.Search_Type == 0);

TA_prev_TP_Final_Eff_Monkey3 = TA_prev_TP_Final(Eff_TP_Monkey3, :);
TA_prev_TP_Final_Ineff_Monkey3 = TA_prev_TP_Final(Ineff_TP_Monkey3, :);

% Monkey 4
Eff_TP_Monkey4 = find(TA_prev_TP_Final.Subject == 4 & TA_prev_TP_Final.Search_Type == 1);
Ineff_TP_Monkey4 = find(TA_prev_TP_Final.Subject == 4 & TA_prev_TP_Final.Search_Type == 0);

TA_prev_TP_Final_Eff_Monkey4 = TA_prev_TP_Final(Eff_TP_Monkey4, :);
TA_prev_TP_Final_Ineff_Monkey4 = TA_prev_TP_Final(Ineff_TP_Monkey4, :);

% Separate Eff and Ineff conditions for TA_prev_TA_Final for each monkey
% Monkey 3
Eff_TA_Monkey3 = find(TA_prev_TA_Final.Subject == 3 & TA_prev_TA_Final.Search_Type == 1);
Ineff_TA_Monkey3 = find(TA_prev_TA_Final.Subject == 3 & TA_prev_TA_Final.Search_Type == 0);

TA_prev_TA_Final_Eff_Monkey3 = TA_prev_TA_Final(Eff_TA_Monkey3, :);
TA_prev_TA_Final_Ineff_Monkey3 = TA_prev_TA_Final(Ineff_TA_Monkey3, :);

% Monkey 4
Eff_TA_Monkey4 = find(TA_prev_TA_Final.Subject == 4 & TA_prev_TA_Final.Search_Type == 1);
Ineff_TA_Monkey4 = find(TA_prev_TA_Final.Subject == 4 & TA_prev_TA_Final.Search_Type == 0);

TA_prev_TA_Final_Eff_Monkey4 = TA_prev_TA_Final(Eff_TA_Monkey4, :);
TA_prev_TA_Final_Ineff_Monkey4 = TA_prev_TA_Final(Ineff_TA_Monkey4, :);

% Display size separation for Eff and Ineff for Monkey 3 and Monkey 4 for TA Prev TP

% TA_TP display size separation for Eff for Monkey 3
DspSize3_TP_Eff_Monkey3 = find(TA_prev_TP_Final_Eff_Monkey3.DispSize == 3);
DspSize5_TP_Eff_Monkey3 = find(TA_prev_TP_Final_Eff_Monkey3.DispSize == 5);
DspSize7_TP_Eff_Monkey3 = find(TA_prev_TP_Final_Eff_Monkey3.DispSize == 7);
DspSize9_TP_Eff_Monkey3 = find(TA_prev_TP_Final_Eff_Monkey3.DispSize == 9);

Monkey3_TP_Eff_SearchTime_3 = nanmean(TA_prev_TP_Final_Eff_Monkey3(DspSize3_TP_Eff_Monkey3, :).SearchTime);
Monkey3_TP_Eff_SearchTime_5 = nanmean(TA_prev_TP_Final_Eff_Monkey3(DspSize5_TP_Eff_Monkey3, :).SearchTime);
Monkey3_TP_Eff_SearchTime_7 = nanmean(TA_prev_TP_Final_Eff_Monkey3(DspSize7_TP_Eff_Monkey3, :).SearchTime);
Monkey3_TP_Eff_SearchTime_9 = nanmean(TA_prev_TP_Final_Eff_Monkey3(DspSize9_TP_Eff_Monkey3, :).SearchTime);

Monkey3_TP_Eff_SEM_3 = nanstd(TA_prev_TP_Final_Eff_Monkey3(DspSize3_TP_Eff_Monkey3, :).SearchTime) / sqrt(length(DspSize3_TP_Eff_Monkey3));
Monkey3_TP_Eff_SEM_5 = nanstd(TA_prev_TP_Final_Eff_Monkey3(DspSize5_TP_Eff_Monkey3, :).SearchTime) / sqrt(length(DspSize5_TP_Eff_Monkey3));
Monkey3_TP_Eff_SEM_7 = nanstd(TA_prev_TP_Final_Eff_Monkey3(DspSize7_TP_Eff_Monkey3, :).SearchTime) / sqrt(length(DspSize7_TP_Eff_Monkey3));
Monkey3_TP_Eff_SEM_9 = nanstd(TA_prev_TP_Final_Eff_Monkey3(DspSize9_TP_Eff_Monkey3, :).SearchTime) / sqrt(length(DspSize9_TP_Eff_Monkey3));

% Repeat the same for Ineff for Monkey 3
DspSize3_TP_Ineff_Monkey3 = find(TA_prev_TP_Final_Ineff_Monkey3.DispSize == 3);
DspSize5_TP_Ineff_Monkey3 = find(TA_prev_TP_Final_Ineff_Monkey3.DispSize == 5);
DspSize7_TP_Ineff_Monkey3 = find(TA_prev_TP_Final_Ineff_Monkey3.DispSize == 7);
DspSize9_TP_Ineff_Monkey3 = find(TA_prev_TP_Final_Ineff_Monkey3.DispSize == 9);

Monkey3_TP_Ineff_SearchTime_3 = nanmean(TA_prev_TP_Final_Ineff_Monkey3(DspSize3_TP_Ineff_Monkey3, :).SearchTime);
Monkey3_TP_Ineff_SearchTime_5 = nanmean(TA_prev_TP_Final_Ineff_Monkey3(DspSize5_TP_Ineff_Monkey3, :).SearchTime);
Monkey3_TP_Ineff_SearchTime_7 = nanmean(TA_prev_TP_Final_Ineff_Monkey3(DspSize7_TP_Ineff_Monkey3, :).SearchTime);
Monkey3_TP_Ineff_SearchTime_9 = nanmean(TA_prev_TP_Final_Ineff_Monkey3(DspSize9_TP_Ineff_Monkey3, :).SearchTime);

Monkey3_TP_Ineff_SEM_3 = nanstd(TA_prev_TP_Final_Ineff_Monkey3(DspSize3_TP_Ineff_Monkey3, :).SearchTime) / sqrt(length(DspSize3_TP_Ineff_Monkey3));
Monkey3_TP_Ineff_SEM_5 = nanstd(TA_prev_TP_Final_Ineff_Monkey3(DspSize5_TP_Ineff_Monkey3, :).SearchTime) / sqrt(length(DspSize5_TP_Ineff_Monkey3));
Monkey3_TP_Ineff_SEM_7 = nanstd(TA_prev_TP_Final_Ineff_Monkey3(DspSize7_TP_Ineff_Monkey3, :).SearchTime) / sqrt(length(DspSize7_TP_Ineff_Monkey3));
Monkey3_TP_Ineff_SEM_9 = nanstd(TA_prev_TP_Final_Ineff_Monkey3(DspSize9_TP_Ineff_Monkey3, :).SearchTime) / sqrt(length(DspSize9_TP_Ineff_Monkey3));

% Display size separation for Eff and Ineff for Monkey 4

% TA_TP display size separation for Eff for Monkey 4
DspSize3_TP_Eff_Monkey4 = find(TA_prev_TP_Final_Eff_Monkey4.DispSize == 3);
DspSize5_TP_Eff_Monkey4 = find(TA_prev_TP_Final_Eff_Monkey4.DispSize == 5);
DspSize7_TP_Eff_Monkey4 = find(TA_prev_TP_Final_Eff_Monkey4.DispSize == 7);
DspSize9_TP_Eff_Monkey4 = find(TA_prev_TP_Final_Eff_Monkey4.DispSize == 9);

Monkey4_TP_Eff_SearchTime_3 = nanmean(TA_prev_TP_Final_Eff_Monkey4(DspSize3_TP_Eff_Monkey4, :).SearchTime);
Monkey4_TP_Eff_SearchTime_5 = nanmean(TA_prev_TP_Final_Eff_Monkey4(DspSize5_TP_Eff_Monkey4, :).SearchTime);
Monkey4_TP_Eff_SearchTime_7 = nanmean(TA_prev_TP_Final_Eff_Monkey4(DspSize7_TP_Eff_Monkey4, :).SearchTime);
Monkey4_TP_Eff_SearchTime_9 = nanmean(TA_prev_TP_Final_Eff_Monkey4(DspSize9_TP_Eff_Monkey4, :).SearchTime);

Monkey4_TP_Eff_SEM_3 = nanstd(TA_prev_TP_Final_Eff_Monkey4(DspSize3_TP_Eff_Monkey4, :).SearchTime) / sqrt(length(DspSize3_TP_Eff_Monkey4));
Monkey4_TP_Eff_SEM_5 = nanstd(TA_prev_TP_Final_Eff_Monkey4(DspSize5_TP_Eff_Monkey4, :).SearchTime) / sqrt(length(DspSize5_TP_Eff_Monkey4));
Monkey4_TP_Eff_SEM_7 = nanstd(TA_prev_TP_Final_Eff_Monkey4(DspSize7_TP_Eff_Monkey4, :).SearchTime) / sqrt(length(DspSize7_TP_Eff_Monkey4));
Monkey4_TP_Eff_SEM_9 = nanstd(TA_prev_TP_Final_Eff_Monkey4(DspSize9_TP_Eff_Monkey4, :).SearchTime) / sqrt(length(DspSize9_TP_Eff_Monkey4));

% Repeat the same for Ineff for Monkey 4
DspSize3_TP_Ineff_Monkey4 = find(TA_prev_TP_Final_Ineff_Monkey4.DispSize == 3);
DspSize5_TP_Ineff_Monkey4 = find(TA_prev_TP_Final_Ineff_Monkey4.DispSize == 5);
DspSize7_TP_Ineff_Monkey4 = find(TA_prev_TP_Final_Ineff_Monkey4.DispSize == 7);
DspSize9_TP_Ineff_Monkey4 = find(TA_prev_TP_Final_Ineff_Monkey4.DispSize == 9);

Monkey4_TP_Ineff_SearchTime_3 = nanmean(TA_prev_TP_Final_Ineff_Monkey4(DspSize3_TP_Ineff_Monkey4, :).SearchTime);
Monkey4_TP_Ineff_SearchTime_5 = nanmean(TA_prev_TP_Final_Ineff_Monkey4(DspSize5_TP_Ineff_Monkey4, :).SearchTime);
Monkey4_TP_Ineff_SearchTime_7 = nanmean(TA_prev_TP_Final_Ineff_Monkey4(DspSize7_TP_Ineff_Monkey4, :).SearchTime);
Monkey4_TP_Ineff_SearchTime_9 = nanmean(TA_prev_TP_Final_Ineff_Monkey4(DspSize9_TP_Ineff_Monkey4, :).SearchTime);

Monkey4_TP_Ineff_SEM_3 = nanstd(TA_prev_TP_Final_Ineff_Monkey4(DspSize3_TP_Ineff_Monkey4, :).SearchTime) / sqrt(length(DspSize3_TP_Ineff_Monkey4));
Monkey4_TP_Ineff_SEM_5 = nanstd(TA_prev_TP_Final_Ineff_Monkey4(DspSize5_TP_Ineff_Monkey4, :).SearchTime) / sqrt(length(DspSize5_TP_Ineff_Monkey4));
Monkey4_TP_Ineff_SEM_7 = nanstd(TA_prev_TP_Final_Ineff_Monkey4(DspSize7_TP_Ineff_Monkey4, :).SearchTime) / sqrt(length(DspSize7_TP_Ineff_Monkey4));
Monkey4_TP_Ineff_SEM_9 = nanstd(TA_prev_TP_Final_Ineff_Monkey4(DspSize9_TP_Ineff_Monkey4, :).SearchTime) / sqrt(length(DspSize9_TP_Ineff_Monkey4));

% Repeat for TA_prev_TA_Final for both monkeys

% TA_TA display size separation for Eff for Monkey 3
DspSize3_TA_Eff_Monkey3 = find(TA_prev_TA_Final_Eff_Monkey3.DispSize == 3);
DspSize5_TA_Eff_Monkey3 = find(TA_prev_TA_Final_Eff_Monkey3.DispSize == 5);
DspSize7_TA_Eff_Monkey3 = find(TA_prev_TA_Final_Eff_Monkey3.DispSize == 7);
DspSize9_TA_Eff_Monkey3 = find(TA_prev_TA_Final_Eff_Monkey3.DispSize == 9);

Monkey3_TA_Eff_SearchTime_3 = nanmean(TA_prev_TA_Final_Eff_Monkey3(DspSize3_TA_Eff_Monkey3, :).SearchTime);
Monkey3_TA_Eff_SearchTime_5 = nanmean(TA_prev_TA_Final_Eff_Monkey3(DspSize5_TA_Eff_Monkey3, :).SearchTime);
Monkey3_TA_Eff_SearchTime_7 = nanmean(TA_prev_TA_Final_Eff_Monkey3(DspSize7_TA_Eff_Monkey3, :).SearchTime);
Monkey3_TA_Eff_SearchTime_9 = nanmean(TA_prev_TA_Final_Eff_Monkey3(DspSize9_TA_Eff_Monkey3, :).SearchTime);

Monkey3_TA_Eff_SEM_3 = nanstd(TA_prev_TA_Final_Eff_Monkey3(DspSize3_TA_Eff_Monkey3, :).SearchTime) / sqrt(length(DspSize3_TA_Eff_Monkey3));
Monkey3_TA_Eff_SEM_5 = nanstd(TA_prev_TA_Final_Eff_Monkey3(DspSize5_TA_Eff_Monkey3, :).SearchTime) / sqrt(length(DspSize5_TA_Eff_Monkey3));
Monkey3_TA_Eff_SEM_7 = nanstd(TA_prev_TA_Final_Eff_Monkey3(DspSize7_TA_Eff_Monkey3, :).SearchTime) / sqrt(length(DspSize7_TA_Eff_Monkey3));
Monkey3_TA_Eff_SEM_9 = nanstd(TA_prev_TA_Final_Eff_Monkey3(DspSize9_TA_Eff_Monkey3, :).SearchTime) / sqrt(length(DspSize9_TA_Eff_Monkey3));

% Repeat the same for Ineff for Monkey 3
DspSize3_TA_Ineff_Monkey3 = find(TA_prev_TA_Final_Ineff_Monkey3.DispSize == 3);
DspSize5_TA_Ineff_Monkey3 = find(TA_prev_TA_Final_Ineff_Monkey3.DispSize == 5);
DspSize7_TA_Ineff_Monkey3 = find(TA_prev_TA_Final_Ineff_Monkey3.DispSize == 7);
DspSize9_TA_Ineff_Monkey3 = find(TA_prev_TA_Final_Ineff_Monkey3.DispSize == 9);

Monkey3_TA_Ineff_SearchTime_3 = nanmean(TA_prev_TA_Final_Ineff_Monkey3(DspSize3_TA_Ineff_Monkey3, :).SearchTime);
Monkey3_TA_Ineff_SearchTime_5 = nanmean(TA_prev_TA_Final_Ineff_Monkey3(DspSize5_TA_Ineff_Monkey3, :).SearchTime);
Monkey3_TA_Ineff_SearchTime_7 = nanmean(TA_prev_TA_Final_Ineff_Monkey3(DspSize7_TA_Ineff_Monkey3, :).SearchTime);
Monkey3_TA_Ineff_SearchTime_9 = nanmean(TA_prev_TA_Final_Ineff_Monkey3(DspSize9_TA_Ineff_Monkey3, :).SearchTime);

Monkey3_TA_Ineff_SEM_3 = nanstd(TA_prev_TA_Final_Ineff_Monkey3(DspSize3_TA_Ineff_Monkey3, :).SearchTime) / sqrt(length(DspSize3_TA_Ineff_Monkey3));
Monkey3_TA_Ineff_SEM_5 = nanstd(TA_prev_TA_Final_Ineff_Monkey3(DspSize5_TA_Ineff_Monkey3, :).SearchTime) / sqrt(length(DspSize5_TA_Ineff_Monkey3));
Monkey3_TA_Ineff_SEM_7 = nanstd(TA_prev_TA_Final_Ineff_Monkey3(DspSize7_TA_Ineff_Monkey3, :).SearchTime) / sqrt(length(DspSize7_TA_Ineff_Monkey3));
Monkey3_TA_Ineff_SEM_9 = nanstd(TA_prev_TA_Final_Ineff_Monkey3(DspSize9_TA_Ineff_Monkey3, :).SearchTime) / sqrt(length(DspSize9_TA_Ineff_Monkey3));

% Display size separation for Eff and Ineff for Monkey 4

% TA_TA display size separation for Eff for Monkey 4
DspSize3_TA_Eff_Monkey4 = find(TA_prev_TA_Final_Eff_Monkey4.DispSize == 3);
DspSize5_TA_Eff_Monkey4 = find(TA_prev_TA_Final_Eff_Monkey4.DispSize == 5);
DspSize7_TA_Eff_Monkey4 = find(TA_prev_TA_Final_Eff_Monkey4.DispSize == 7);
DspSize9_TA_Eff_Monkey4 = find(TA_prev_TA_Final_Eff_Monkey4.DispSize == 9);

Monkey4_TA_Eff_SearchTime_3 = nanmean(TA_prev_TA_Final_Eff_Monkey4(DspSize3_TA_Eff_Monkey4, :).SearchTime);
Monkey4_TA_Eff_SearchTime_5 = nanmean(TA_prev_TA_Final_Eff_Monkey4(DspSize5_TA_Eff_Monkey4, :).SearchTime);
Monkey4_TA_Eff_SearchTime_7 = nanmean(TA_prev_TA_Final_Eff_Monkey4(DspSize7_TA_Eff_Monkey4, :).SearchTime);
Monkey4_TA_Eff_SearchTime_9 = nanmean(TA_prev_TA_Final_Eff_Monkey4(DspSize9_TA_Eff_Monkey4, :).SearchTime);

Monkey4_TA_Eff_SEM_3 = nanstd(TA_prev_TA_Final_Eff_Monkey4(DspSize3_TA_Eff_Monkey4, :).SearchTime) / sqrt(length(DspSize3_TA_Eff_Monkey4));
Monkey4_TA_Eff_SEM_5 = nanstd(TA_prev_TA_Final_Eff_Monkey4(DspSize5_TA_Eff_Monkey4, :).SearchTime) / sqrt(length(DspSize5_TA_Eff_Monkey4));
Monkey4_TA_Eff_SEM_7 = nanstd(TA_prev_TA_Final_Eff_Monkey4(DspSize7_TA_Eff_Monkey4, :).SearchTime) / sqrt(length(DspSize7_TA_Eff_Monkey4));
Monkey4_TA_Eff_SEM_9 = nanstd(TA_prev_TA_Final_Eff_Monkey4(DspSize9_TA_Eff_Monkey4, :).SearchTime) / sqrt(length(DspSize9_TA_Eff_Monkey4));

% Repeat the same for Ineff for Monkey 4
DspSize3_TA_Ineff_Monkey4 = find(TA_prev_TA_Final_Ineff_Monkey4.DispSize == 3);
DspSize5_TA_Ineff_Monkey4 = find(TA_prev_TA_Final_Ineff_Monkey4.DispSize == 5);
DspSize7_TA_Ineff_Monkey4 = find(TA_prev_TA_Final_Ineff_Monkey4.DispSize == 7);
DspSize9_TA_Ineff_Monkey4 = find(TA_prev_TA_Final_Ineff_Monkey4.DispSize == 9);

Monkey4_TA_Ineff_SearchTime_3 = nanmean(TA_prev_TA_Final_Ineff_Monkey4(DspSize3_TA_Ineff_Monkey4, :).SearchTime);
Monkey4_TA_Ineff_SearchTime_5 = nanmean(TA_prev_TA_Final_Ineff_Monkey4(DspSize5_TA_Ineff_Monkey4, :).SearchTime);
Monkey4_TA_Ineff_SearchTime_7 = nanmean(TA_prev_TA_Final_Ineff_Monkey4(DspSize7_TA_Ineff_Monkey4, :).SearchTime);
Monkey4_TA_Ineff_SearchTime_9 = nanmean(TA_prev_TA_Final_Ineff_Monkey4(DspSize9_TA_Ineff_Monkey4, :).SearchTime);

Monkey4_TA_Ineff_SEM_3 = nanstd(TA_prev_TA_Final_Ineff_Monkey4(DspSize3_TA_Ineff_Monkey4, :).SearchTime) / sqrt(length(DspSize3_TA_Ineff_Monkey4));
Monkey4_TA_Ineff_SEM_5 = nanstd(TA_prev_TA_Final_Ineff_Monkey4(DspSize5_TA_Ineff_Monkey4, :).SearchTime) / sqrt(length(DspSize5_TA_Ineff_Monkey4));
Monkey4_TA_Ineff_SEM_7 = nanstd(TA_prev_TA_Final_Ineff_Monkey4(DspSize7_TA_Ineff_Monkey4, :).SearchTime) / sqrt(length(DspSize7_TA_Ineff_Monkey4));
Monkey4_TA_Ineff_SEM_9 = nanstd(TA_prev_TA_Final_Ineff_Monkey4(DspSize9_TA_Ineff_Monkey4, :).SearchTime) / sqrt(length(DspSize9_TA_Ineff_Monkey4));

% Prepare data for plotting
display_sizes = [3, 5, 7, 9];

% Prepare data for plotting
display_sizes = [3, 5, 7, 9];

% Plot with error bars for Ineff
figure;

% TP Prev TP Ineff 
errorbar(display_sizes, [Monkey3_TP_Ineff_SearchTime_3, Monkey3_TP_Ineff_SearchTime_5, Monkey3_TP_Ineff_SearchTime_7, Monkey3_TP_Ineff_SearchTime_9], ...
    [Monkey3_TP_Ineff_SEM_3, Monkey3_TP_Ineff_SEM_5, Monkey3_TP_Ineff_SEM_7, Monkey3_TP_Ineff_SEM_9], 'o-', 'LineWidth', 2, 'DisplayName', 'TA Prev TP Ineff', 'Color', 'r');
hold on;
errorbar(display_sizes, [Monkey3_TA_Ineff_SearchTime_3, Monkey3_TA_Ineff_SearchTime_5, Monkey3_TA_Ineff_SearchTime_7, Monkey3_TA_Ineff_SearchTime_9], ...
    [Monkey3_TA_Ineff_SEM_3, Monkey3_TA_Ineff_SEM_5, Monkey3_TA_Ineff_SEM_7, Monkey3_TA_Ineff_SEM_9], 'o-', 'LineWidth', 2, 'DisplayName', 'TA Prev TA Ineff', 'Color', 'b');
xlabel('Display Size');
ylabel('Mean Search Time');
title('Mean Search Time by Display Size for Ineff (Monkeys 3)');
legend('show');
grid on;
hold off;

% Repeat similar plotting for Eff conditions
figure;

% TP Prev TP Eff (Monkey 3 in red, Monkey 4 in blue)
errorbar(display_sizes, [Monkey3_TP_Eff_SearchTime_3, Monkey3_TP_Eff_SearchTime_5, Monkey3_TP_Eff_SearchTime_7, Monkey3_TP_Eff_SearchTime_9], ...
    [Monkey3_TP_Eff_SEM_3, Monkey3_TP_Eff_SEM_5, Monkey3_TP_Eff_SEM_7, Monkey3_TP_Eff_SEM_9], 'o-', 'LineWidth', 2, 'DisplayName', 'TA Prev TP Eff', 'Color', 'r');
hold on;
errorbar(display_sizes, [Monkey3_TA_Eff_SearchTime_3, Monkey3_TA_Eff_SearchTime_5, Monkey3_TA_Eff_SearchTime_7, Monkey3_TA_Eff_SearchTime_9], ...
    [Monkey3_TA_Eff_SEM_3, Monkey3_TA_Eff_SEM_5, Monkey3_TA_Eff_SEM_7, Monkey3_TA_Eff_SEM_9], 'o-', 'LineWidth', 2, 'DisplayName', 'TA Prev TA Eff', 'Color', 'b');
xlabel('Display Size');
ylabel('Mean Search Time');
title('Mean Search Time by Display Size for Eff (Monkeys 3)');
legend('show');
grid on;
hold off;

% Repeat similar plotting for TP Prev TA Eff and Ineff conditions
figure;

% TP Prev TA Ineff (Monkey 3 in red, Monkey 4 in blue)
errorbar(display_sizes, [Monkey4_TP_Ineff_SearchTime_3, Monkey4_TP_Ineff_SearchTime_5, Monkey4_TP_Ineff_SearchTime_7, Monkey4_TP_Ineff_SearchTime_9], ...
    [Monkey4_TP_Ineff_SEM_3, Monkey4_TP_Ineff_SEM_5, Monkey4_TP_Ineff_SEM_7, Monkey4_TP_Ineff_SEM_9], 's-', 'LineWidth', 2, 'DisplayName', 'TA Prev TP Ineff', 'Color', 'r');
hold on;
errorbar(display_sizes, [Monkey4_TA_Ineff_SearchTime_3, Monkey4_TA_Ineff_SearchTime_5, Monkey4_TA_Ineff_SearchTime_7, Monkey4_TA_Ineff_SearchTime_9], ...
    [Monkey4_TA_Ineff_SEM_3, Monkey4_TA_Ineff_SEM_5, Monkey4_TA_Ineff_SEM_7, Monkey4_TA_Ineff_SEM_9], 's-', 'LineWidth', 2, 'DisplayName', 'TA Prev TA Ineff', 'Color', 'b');
xlabel('Display Size');
ylabel('Mean Search Time');
title('Mean Search Time by Display Size for Ineff (Monkey 4)');
legend('show');
grid on;
hold off;

figure;

% TP Prev TA Eff (Monkey 3 in red, Monkey 4 in blue)
errorbar(display_sizes, [Monkey4_TP_Eff_SearchTime_3, Monkey4_TP_Eff_SearchTime_5, Monkey4_TP_Eff_SearchTime_7, Monkey4_TP_Eff_SearchTime_9], ...
    [Monkey4_TP_Eff_SEM_3, Monkey4_TP_Eff_SEM_5, Monkey4_TP_Eff_SEM_7, Monkey4_TP_Eff_SEM_9], 's-', 'LineWidth', 2, 'DisplayName', 'TA Prev TP Eff', 'Color', 'r');
hold on;
errorbar(display_sizes, [Monkey4_TA_Eff_SearchTime_3, Monkey4_TA_Eff_SearchTime_5, Monkey4_TA_Eff_SearchTime_7, Monkey4_TA_Eff_SearchTime_9], ...
    [Monkey4_TA_Eff_SEM_3, Monkey4_TA_Eff_SEM_5, Monkey4_TA_Eff_SEM_7, Monkey4_TA_Eff_SEM_9], 's-', 'LineWidth', 2, 'DisplayName', 'TA Prev TA Eff', 'Color', 'b');
xlabel('Display Size');
ylabel('Mean Search Time');
title('Mean Search Time by Display Size for Eff (Monkey 4)');
legend('show');
grid on;
hold off;








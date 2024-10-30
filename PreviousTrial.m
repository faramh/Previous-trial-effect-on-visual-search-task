% Is there any effect for previous trial on the current trial?
% Lets answer this Question with LFP analysis

%% Load The data for SNr
clc; clear; close all;

%load('complete_data_both.mat');

% vlPFC
%load('data_both_final_tptfnewmh.mat');

% SNr
load('data_JP_final_tptf_12_18.mat');
%% Condition separation
% C= Current P=Previous
% Conditions: C:TP/P:TP -- C:TP/P:TA -- C:TA/P:TP -- C:TA/P:TA


% Get all unique session names
sessions = unique(data.Session_Name);

% Initialize arrays to hold the concatenated results
TP_prev_TP = [];
TP_prev_TA = [];
TA_prev_TP = [];
TA_prev_TA = [];

TP_prev_TP_Pre = [];
TP_prev_TA_Pre = [];
TA_prev_TP_Pre = [];
TA_prev_TA_Pre = [];

counter=0;
% Loop through each session
for i = 1:length(sessions)
    % Extract the data for the current session
    current_session = sessions{i};
    session_data = data(strcmp(data.Session_Name, current_session), :);
    
    % Sort the session data by Trial_N to ensure sequential order
    % [~, idx] = sort(session_data.Trial_N);
    % session_data = session_data(idx, :);

    TrialsIndices=session_data.Trial_N;
    
    % Loop through trials, starting from the second one (n-1 concept)
    for j = 2:height(session_data)

        current_trial = session_data(j, :);
        previous_trial = session_data(j-1, :);

        CurrentTrialNumber=current_trial.Trial_N;
        PrevTrialNumber=CurrentTrialNumber-1;

        %counter=counter+1;

        
        if ismember(PrevTrialNumber, TrialsIndices)
            %disp(['Current: ', num2str(current_trial.Trial_N), ' - Previous: ', num2str(previous_trial.Trial_N)]);
        % Check current and previous trial types
            if current_trial.TA_TP == 4 && previous_trial.TA_TP == 4
                % TP and previous TP
                %disp(['Current is:',num2str(current_trial.Trial_N), '  And prev is ',num2str(previous_trial.Trial_N)]);
                TP_prev_TP = [TP_prev_TP; current_trial];
                TP_prev_TP_Pre=[TP_prev_TP_Pre; previous_trial];
            elseif current_trial.TA_TP == 4 && previous_trial.TA_TP == 3
                % TP and previous TA
                TP_prev_TA = [TP_prev_TA; current_trial];
                TP_prev_TA_Pre=[TP_prev_TA_Pre; previous_trial];
                
            elseif current_trial.TA_TP == 3 && previous_trial.TA_TP == 4
                % TA and previous TP
                TA_prev_TP = [TA_prev_TP; current_trial];
                TA_prev_TP_Pre=[TA_prev_TP_Pre; previous_trial];
                
            elseif current_trial.TA_TP == 3 && previous_trial.TA_TP == 3
                % TA and previous TA
                TA_prev_TA = [TA_prev_TA; current_trial];
                TA_prev_TA_Pre=[TA_prev_TA_Pre; previous_trial];
                
            end
        else
            counter=counter+1;
           
        end
    end
end

% Now TP_prev_TP, TP_prev_TA, TA_prev_TP, TA_prev_TA contain the trials
% from all sessions concatenated in the respective categories.




%% Extracting LFP DATA from Table
wname = 'amor';
t_baseline = 0.5;
Fs=1000;

time = -0.5:1/Fs:3-1/Fs;

TP_prev_TP_Data=cell2mat((TP_prev_TP.Data)');

TP_prev_TA_Data=cell2mat((TP_prev_TA.Data)');

TA_prev_TP_Data=cell2mat((TA_prev_TP.Data)');

TA_prev_TA_Data=cell2mat((TA_prev_TA.Data)');

[wt_TP_prev_TA_avg, f_wt_TP_prev_TA_avg] = cal_cwt_1(TP_prev_TA_Data', wname, Fs,3, 'TP eff Corr');
wt_TP_prev_TA_avg = baseline_normalization_mat(wt_TP_prev_TA_avg, t_baseline, Fs);

[wt_TA_prev_TA_avg, f_wt_TA_prev_TA_avg] = cal_cwt_1(TA_prev_TA_Data', wname, Fs,3, 'TP eff Corr');
wt_TA_prev_TA_avg = baseline_normalization_mat(wt_TA_prev_TA_avg, t_baseline, Fs);

[wt_TA_prev_TP_avg, f_wt_TA_prev_TP_avg] = cal_cwt_1(TA_prev_TP_Data', wname, Fs,3, 'TP eff Corr');
wt_TA_prev_TP_avg = baseline_normalization_mat(wt_TA_prev_TP_avg, t_baseline, Fs);

[wt_TP_prev_TP_avg, f_wt_TP_prev_TP_avg] = cal_cwt_1(TP_prev_TP_Data', wname, Fs,3, 'TP eff Corr');
wt_TP_prev_TP_avg = baseline_normalization_mat(wt_TP_prev_TP_avg, t_baseline, Fs);


%% Smooth matrices across time

windowSize=35;
wt_TP_prev_TA_avg_smoothed = smoothdata(wt_TP_prev_TA_avg, 2, 'gaussian', windowSize); % Smooth along the time axis (dimension 2)

wt_TP_prev_TP_avg_smoothed = smoothdata(wt_TP_prev_TP_avg, 2, 'gaussian', windowSize); % Smooth along the time axis (dimension 2)

wt_TA_prev_TA_avg_smoothed = smoothdata(wt_TA_prev_TA_avg, 2, 'gaussian', windowSize); % Smooth along the time axis (dimension 2)

wt_TA_prev_TP_avg_smoothed = smoothdata(wt_TA_prev_TP_avg, 2, 'gaussian', windowSize); % Smooth along the time axis (dimension 2)


%% Plot LFPs

fig = figure('units','normalized','outerposition',[0 0 1 1]);

subplot(2,2,1);
pcolor(time,f_wt_TP_prev_TP_avg,wt_TP_prev_TP_avg_smoothed);
%xline(0,'LineWidth',1.5,'LineStyle','--'); % baseline
set(gca,'FontSize',15); 

title("Current:Target Present - Previous: Target Present")
xlabel('time(s)','interpreter','latex');
ylabel('frequency(Hz)','interpreter','latex');
shading flat;
ylim([0 200]);
xlim([-0.5,2.5]);
colormap jet
caxis([-7 7]);


subplot(2,2,2);

pcolor(time,f_wt_TP_prev_TA_avg,wt_TP_prev_TA_avg_smoothed);
%xline(0,'LineWidth',1.5,'LineStyle','--'); % baseline
set(gca,'FontSize',15); 

title("Current:Target Present - Previous: Target Absent")
xlabel('time(s)','interpreter','latex');
ylabel('frequency(Hz)','interpreter','latex');
shading flat;
ylim([0 200]);
xlim([-0.5,2.5]);
colormap jet
caxis([-7 7]);
 
subplot(2,2,3);

pcolor(time,f_wt_TA_prev_TP_avg,wt_TA_prev_TP_avg_smoothed);
%xline(0,'LineWidth',1.5,'LineStyle','--'); % baseline
set(gca,'FontSize',15); 

title("Current:Target Absent - Previous: Target Present")
xlabel('time(s)','interpreter','latex');
ylabel('frequency(Hz)','interpreter','latex');
shading flat;
ylim([0 200]);
xlim([-0.5,2.5]);
colormap jet
caxis([-7 7]);


subplot(2,2,4);

pcolor(time,f_wt_TA_prev_TA_avg,wt_TA_prev_TA_avg_smoothed);
%xline(0,'LineWidth',1.5,'LineStyle','--'); % baseline
set(gca,'FontSize',15); 

title("Current:Target Absent - Previous: Target Absent")
xlabel('time(s)','interpreter','latex');
ylabel('frequency(Hz)','interpreter','latex');
shading flat;
ylim([0 200]);
xlim([-0.5,2.5]);
colormap jet
caxis([-7 7]);


%% Considering EFF/INEFF conditions 


% TA - TA
TA_prev_TA_EFF_temp = find(TA_prev_TA.Search_Type=="p");
TA_prev_TA_InEFF_temp = find(TA_prev_TA.Search_Type=="s");

TA_prev_TA_EFF=TA_prev_TA(TA_prev_TA_EFF_temp,:);
TA_prev_TA_InEFF=TA_prev_TA(TA_prev_TA_InEFF_temp,:);

% TA - TP
TA_prev_TP_EFF_temp = find(TA_prev_TP.Search_Type=="p");
TA_prev_TP_InEFF_temp = find(TA_prev_TP.Search_Type=="s");

TA_prev_TP_EFF=TA_prev_TP(TA_prev_TP_EFF_temp,:);
TA_prev_TP_InEFF=TA_prev_TP(TA_prev_TP_InEFF_temp,:);


% TP - TP
TP_prev_TP_EFF_temp = find(TP_prev_TP.Search_Type=="p");
TP_prev_TP_InEFF_temp = find(TP_prev_TP.Search_Type=="s");

TP_prev_TP_EFF=TP_prev_TP(TP_prev_TP_EFF_temp,:);
TP_prev_TP_InEFF=TP_prev_TP(TP_prev_TP_InEFF_temp,:);



% TP - TA
TP_prev_TA_EFF_temp = find(TP_prev_TA.Search_Type=="p");
TP_prev_TA_InEFF_temp = find(TP_prev_TA.Search_Type=="s");

TP_prev_TA_EFF=TP_prev_TA(TP_prev_TA_EFF_temp,:);
TP_prev_TA_InEFF=TP_prev_TA(TP_prev_TA_InEFF_temp,:);



%% Extracting LFP DATA from Table
wname = 'amor';
t_baseline = 0.5;
Fs=1000;

time = -0.5:1/Fs:3-1/Fs;

% Convert cell arrays to matrices for EFF and InEFF sections
TP_prev_TP_Data_EFF = cell2mat((TP_prev_TP_EFF.Data)');
TP_prev_TP_Data_InEFF = cell2mat((TP_prev_TP_InEFF.Data)');

TP_prev_TA_Data_EFF = cell2mat((TP_prev_TA_EFF.Data)');
TP_prev_TA_Data_InEFF = cell2mat((TP_prev_TA_InEFF.Data)');

TA_prev_TP_Data_EFF = cell2mat((TA_prev_TP_EFF.Data)');
TA_prev_TP_Data_InEFF = cell2mat((TA_prev_TP_InEFF.Data)');

TA_prev_TA_Data_EFF = cell2mat((TA_prev_TA_EFF.Data)');
TA_prev_TA_Data_InEFF = cell2mat((TA_prev_TA_InEFF.Data)');

% Calculate wavelet transform and apply baseline normalization for EFF section
[wt_TP_prev_TA_avg_EFF, f_wt_TP_prev_TA_avg_EFF] = cal_cwt_1(TP_prev_TA_Data_EFF', wname, Fs, 3, 'TP eff Corr');
wt_TP_prev_TA_avg_EFF = baseline_normalization_mat(wt_TP_prev_TA_avg_EFF, t_baseline, Fs);

[wt_TA_prev_TA_avg_EFF, f_wt_TA_prev_TA_avg_EFF] = cal_cwt_1(TA_prev_TA_Data_EFF', wname, Fs, 3, 'TP eff Corr');
wt_TA_prev_TA_avg_EFF = baseline_normalization_mat(wt_TA_prev_TA_avg_EFF, t_baseline, Fs);

[wt_TA_prev_TP_avg_EFF, f_wt_TA_prev_TP_avg_EFF] = cal_cwt_1(TA_prev_TP_Data_EFF', wname, Fs, 3, 'TP eff Corr');
wt_TA_prev_TP_avg_EFF = baseline_normalization_mat(wt_TA_prev_TP_avg_EFF, t_baseline, Fs);

[wt_TP_prev_TP_avg_EFF, f_wt_TP_prev_TP_avg_EFF] = cal_cwt_1(TP_prev_TP_Data_EFF', wname, Fs, 3, 'TP eff Corr');
wt_TP_prev_TP_avg_EFF = baseline_normalization_mat(wt_TP_prev_TP_avg_EFF, t_baseline, Fs);

% Calculate wavelet transform and apply baseline normalization for InEFF section
[wt_TP_prev_TA_avg_InEFF, f_wt_TP_prev_TA_avg_InEFF] = cal_cwt_1(TP_prev_TA_Data_InEFF', wname, Fs, 3, 'TP inEff Corr');
wt_TP_prev_TA_avg_InEFF = baseline_normalization_mat(wt_TP_prev_TA_avg_InEFF, t_baseline, Fs);

[wt_TA_prev_TA_avg_InEFF, f_wt_TA_prev_TA_avg_InEFF] = cal_cwt_1(TA_prev_TA_Data_InEFF', wname, Fs, 3, 'TP inEff Corr');
wt_TA_prev_TA_avg_InEFF = baseline_normalization_mat(wt_TA_prev_TA_avg_InEFF, t_baseline, Fs);

[wt_TA_prev_TP_avg_InEFF, f_wt_TA_prev_TP_avg_InEFF] = cal_cwt_1(TA_prev_TP_Data_InEFF', wname, Fs, 3, 'TP inEff Corr');
wt_TA_prev_TP_avg_InEFF = baseline_normalization_mat(wt_TA_prev_TP_avg_InEFF, t_baseline, Fs);

[wt_TP_prev_TP_avg_InEFF, f_wt_TP_prev_TP_avg_InEFF] = cal_cwt_1(TP_prev_TP_Data_InEFF', wname, Fs, 3, 'TP inEff Corr');
wt_TP_prev_TP_avg_InEFF = baseline_normalization_mat(wt_TP_prev_TP_avg_InEFF, t_baseline, Fs);

%% Smooth matrices across time for EFF section
windowSize = 35;
wt_TP_prev_TA_avg_EFF_smoothed = smoothdata(wt_TP_prev_TA_avg_EFF, 2, 'gaussian', windowSize);
wt_TP_prev_TP_avg_EFF_smoothed = smoothdata(wt_TP_prev_TP_avg_EFF, 2, 'gaussian', windowSize);
wt_TA_prev_TA_avg_EFF_smoothed = smoothdata(wt_TA_prev_TA_avg_EFF, 2, 'gaussian', windowSize);
wt_TA_prev_TP_avg_EFF_smoothed = smoothdata(wt_TA_prev_TP_avg_EFF, 2, 'gaussian', windowSize);

%% Smooth matrices across time for InEFF section
wt_TP_prev_TA_avg_InEFF_smoothed = smoothdata(wt_TP_prev_TA_avg_InEFF, 2, 'gaussian', windowSize);
wt_TP_prev_TP_avg_InEFF_smoothed = smoothdata(wt_TP_prev_TP_avg_InEFF, 2, 'gaussian', windowSize);
wt_TA_prev_TA_avg_InEFF_smoothed = smoothdata(wt_TA_prev_TA_avg_InEFF, 2, 'gaussian', windowSize);
wt_TA_prev_TP_avg_InEFF_smoothed = smoothdata(wt_TA_prev_TP_avg_InEFF, 2, 'gaussian', windowSize);



%% Plot considering EFF/InEFF in 2 figures

fig = figure('units','normalized','outerposition',[0 0 1 1]);

% Plotting for EFF Section
subplot(2,2,1);
pcolor(time,f_wt_TP_prev_TP_avg_EFF,wt_TP_prev_TP_avg_EFF_smoothed);
%xline(0,'LineWidth',1.5,'LineStyle','--'); % baseline
set(gca,'FontSize',15); 

title("EFF: Current:Target Present - Previous: Target Present")
xlabel('time(s)','interpreter','latex');
ylabel('frequency(Hz)','interpreter','latex');
shading flat;
ylim([0 200]);
xlim([-0.5,2.5]);
colormap jet
caxis([-7 7]);

subplot(2,2,2);
pcolor(time,f_wt_TP_prev_TA_avg_EFF,wt_TP_prev_TA_avg_EFF_smoothed);
%xline(0,'LineWidth',1.5,'LineStyle','--'); % baseline
set(gca,'FontSize',15); 

title("EFF: Current:Target Present - Previous: Target Absent")
xlabel('time(s)','interpreter','latex');
ylabel('frequency(Hz)','interpreter','latex');
shading flat;
ylim([0 200]);
xlim([-0.5,2.5]);
colormap jet
caxis([-7 7]);

subplot(2,2,3);
pcolor(time,f_wt_TA_prev_TP_avg_EFF,wt_TA_prev_TP_avg_EFF_smoothed);
%xline(0,'LineWidth',1.5,'LineStyle','--'); % baseline
set(gca,'FontSize',15); 

title("EFF: Current:Target Absent - Previous: Target Present")
xlabel('time(s)','interpreter','latex');
ylabel('frequency(Hz)','interpreter','latex');
shading flat;
ylim([0 200]);
xlim([-0.5,2.5]);
colormap jet
caxis([-7 7]);

subplot(2,2,4);
pcolor(time,f_wt_TA_prev_TA_avg_EFF,wt_TA_prev_TA_avg_EFF_smoothed);
%xline(0,'LineWidth',1.5,'LineStyle','--'); % baseline
set(gca,'FontSize',15); 

title("EFF: Current:Target Absent - Previous: Target Absent")
xlabel('time(s)','interpreter','latex');
ylabel('frequency(Hz)','interpreter','latex');
shading flat;
ylim([0 200]);
xlim([-0.5,2.5]);
colormap jet
caxis([-7 7]);

% Create a new figure for InEFF section
fig2 = figure('units','normalized','outerposition',[0 0 1 1]);

% Plotting for InEFF Section
subplot(2,2,1);
pcolor(time,f_wt_TP_prev_TP_avg_InEFF,wt_TP_prev_TP_avg_InEFF_smoothed);
%xline(0,'LineWidth',1.5,'LineStyle','--'); % baseline
set(gca,'FontSize',15); 

title("InEFF: Current:Target Present - Previous: Target Present")
xlabel('time(s)','interpreter','latex');
ylabel('frequency(Hz)','interpreter','latex');
shading flat;
ylim([0 200]);
xlim([-0.5,2.5]);
colormap jet
caxis([-7 7]);

subplot(2,2,2);
pcolor(time,f_wt_TP_prev_TA_avg_InEFF,wt_TP_prev_TA_avg_InEFF_smoothed);
%xline(0,'LineWidth',1.5,'LineStyle','--'); % baseline
set(gca,'FontSize',15); 

title("InEFF: Current:Target Present - Previous: Target Absent")
xlabel('time(s)','interpreter','latex');
ylabel('frequency(Hz)','interpreter','latex');
shading flat;
ylim([0 200]);
xlim([-0.5,2.5]);
colormap jet
caxis([-7 7]);

subplot(2,2,3);
pcolor(time,f_wt_TA_prev_TP_avg_InEFF,wt_TA_prev_TP_avg_InEFF_smoothed);
%xline(0,'LineWidth',1.5,'LineStyle','--'); % baseline
set(gca,'FontSize',15); 

title("InEFF: Current:Target Absent - Previous: Target Present")
xlabel('time(s)','interpreter','latex');
ylabel('frequency(Hz)','interpreter','latex');
shading flat;
ylim([0 200]);
xlim([-0.5,2.5]);
colormap jet
caxis([-7 7]);

subplot(2,2,4);
pcolor(time,f_wt_TA_prev_TA_avg_InEFF,wt_TA_prev_TA_avg_InEFF_smoothed);
%xline(0,'LineWidth',1.5,'LineStyle','--'); % baseline
set(gca,'FontSize',15); 

title("InEFF: Current:Target Absent - Previous: Target Absent")
xlabel('time(s)','interpreter','latex');
ylabel('frequency(Hz)','interpreter','latex');
shading flat;
ylim([0 200]);
xlim([-0.5,2.5]);
colormap jet
caxis([-7 7]);












function [glm_output, encoding_flag, encoding_beta] = stopping_glm_func(neuron_i, mcc_map_info, dataFiles_beh, dirs)

neuralFilename = mcc_map_info.session{neuron_i};
neuronLabel = mcc_map_info.unit{neuron_i};
fprintf('Extracting: %s - %s ... [%i of %i]  \n',neuralFilename,neuronLabel,neuron_i,size(mcc_map_info,1))

%... and find the corresponding behavior file index
behFilename = data_findBehFile(neuralFilename);
behaviorIdx = find(strcmp(dataFiles_beh,behFilename(1:end-4)));

% Load in pre-processed spike data
import_data = struct(); import_data = load_behFile(dirs,behFilename);


%% Extract: get relevant data for GLM table
reg_tbl = table;

for trial_i = 1:size(import_data.events.stateFlags_,1)
    
    % Trial type ---------------------------
    % - No-stop trials
    if import_data.events.stateFlags_.IsGoCorrect(trial_i) == 1
        reg_tbl.trial_type(trial_i) = 1;
    % - Canceled trials
    elseif import_data.events.stateFlags_.IsCancel(trial_i) == 1
        reg_tbl.trial_type(trial_i) = 2;
    % - Non-canceled trials
    elseif import_data.events.stateFlags_.IsNonCancelledNoBrk(trial_i) == 1 ||...
            import_data.events.stateFlags_.IsNonCancelledBrk(trial_i) == 1
        reg_tbl.trial_type(trial_i) = 3;
    % - Abort trials
    else
        reg_tbl.trial_type(trial_i) = 0;
    end
    
    % Stop signal delay ---------------------------
    reg_tbl.ssd(trial_i) = import_data.events.stateFlags_.UseSsdIdx(trial_i);
    
    % Value ---------------------------------------
    if import_data.events.stateFlags_.IsLoRwrd(trial_i) == 1
        reg_tbl.value(trial_i) = 1;
    else
        reg_tbl.value(trial_i) = 2;
    end
    
    % Trial Number ---------------------------------------
    reg_tbl.trial_number(trial_i) = import_data.events.stateFlags_.TrialNumber(trial_i);
    
end

%% Setup spike data into GLM
spk_data_in = load(fullfile('D:\projects\2022-mcc-action\data\','SDF',...
    [neuralFilename '_SDF_' neuronLabel '.mat']));

event_alignment = 'stopSignal_artifical';
baseline_win = [-600:-100];

z_sdf = zscore_sdf(spk_data_in.SDF,[-600:-100],event_alignment);

window_size = 100;
window_shift = 10;
[window_sdf, window_time] = movaverage_sdf(z_sdf, window_size, window_shift);
% 2023-08-26, 23h25: I tested this and it does the job correctly. The
% z-scored profile looks exactly like the raw profile, and the movaverage
% sdf looks like a smoothed version of the raw/zscored sdf.

[n_trials, n_times] = size(window_sdf); % Get the number of windows in the time averaged data

analysis_win_idx = find(window_time >= 0 & window_time <= 600); % Find the relevant indicies for the timepoints of interest
win_fr = nanmean(window_sdf(:,analysis_win_idx),2);

reg_tbl.win_fr = win_fr; % Add the firing rate over this whole window to the GLM table

%% Run GLM: trial type
glm_output.trial_type.sig_times = []; 
glm_output.trial_type.beta_weights = []; 
u_t_mdl = [];

reg_tbl_trialtype = reg_tbl(reg_tbl.trial_type == 1 | reg_tbl.trial_type == 2,:);
   
% For each averaged time point
for timepoint_i = 1:n_times
    
    % Input the timepoint specific firing times
    reg_tbl_trialtype.firing_rate = window_sdf(reg_tbl.trial_type == 1 | reg_tbl.trial_type == 2,timepoint_i);
    
    % Fit the GLM for this
    u_t_mdl = fitlm(reg_tbl_trialtype,'firing_rate ~ trial_type + trial_number');
    
    % GLM output ---------------------------------
    % - Trial type
    glm_output.trial_type.sig_times(1,timepoint_i) = u_t_mdl.Coefficients.pValue(2) < .01; % trial type
    glm_output.trial_type.beta_weights(1,timepoint_i) = u_t_mdl.Coefficients.tStat(2); % trial type
    
    % - Trial number
    glm_output.trial_type.sig_times(2,timepoint_i) = u_t_mdl.Coefficients.pValue(3) < .01; % trial_number
    glm_output.trial_type.beta_weights(2,timepoint_i) = u_t_mdl.Coefficients.tStat(3); % trial_number
    
end

%% Run GLM: stop_signal delay
glm_output.ssd.sig_times = []; 
glm_output.ssd.beta_weights = []; 
u_t_mdl = [];

reg_tbl_stopping = reg_tbl(reg_tbl.trial_type == 2,:);

% For each averaged time point
for timepoint_i = 1:n_times
    
    % Input the timepoint specific firing times
    reg_tbl_stopping.firing_rate = window_sdf(reg_tbl.trial_type == 2,timepoint_i);
    
    % Fit the GLM for this
    u_t_mdl = fitlm(reg_tbl_stopping,'firing_rate ~ ssd*value + trial_number');
    
    % GLM output ---------------------------------   
    % - Stop-signal delay
    glm_output.ssd.sig_times(1,timepoint_i) = u_t_mdl.Coefficients.pValue(2) < .01; % ssd
    glm_output.ssd.beta_weights(1,timepoint_i) = u_t_mdl.Coefficients.tStat(2); % ssd
    
    % - Value
    glm_output.ssd.sig_times(2,timepoint_i) = u_t_mdl.Coefficients.pValue(3) < .01; % value
    glm_output.ssd.beta_weights(2,timepoint_i) = u_t_mdl.Coefficients.tStat(3); % value
        
    % - SSD x Value
    glm_output.ssd.sig_times(3,timepoint_i) = u_t_mdl.Coefficients.pValue(5) < .01; % ssd x value
    glm_output.ssd.beta_weights(3,timepoint_i) = u_t_mdl.Coefficients.tStat(5); % % ssd x value
    
    % - Trial number
    glm_output.ssd.sig_times(4,timepoint_i) = u_t_mdl.Coefficients.pValue(4) < .01; % trial_number
    glm_output.ssd.beta_weights(4,timepoint_i) = u_t_mdl.Coefficients.tStat(4); % trial_number

    
end



%% Determine periods of significance

signal_detect_length = 50;
signal_detect_wins = signal_detect_length/window_shift;

% now ask whether this unit was significant with significance defined as at least
% 100ms of selectivity following SSD

% Trial-type
[~, canc_sig_len, ~] = ZeroOnesCount(glm_output.trial_type.sig_times(1,analysis_win_idx)); % choice direction

% Stop-signal delay
[~, ssd_sig_len, ~] = ZeroOnesCount(glm_output.ssd.sig_times(1,analysis_win_idx)); % choice direction

% Value
[~, value_sig_len, ~] = ZeroOnesCount(glm_output.ssd.sig_times(2,analysis_win_idx)); % choice direction


encoding_flag = []; encoding_beta = [];


encoding_flag(1,1) = any(canc_sig_len >= signal_detect_wins);
encoding_flag(1,2) = any(ssd_sig_len >= signal_detect_wins);
encoding_flag(1,3) = any(value_sig_len >= signal_detect_wins);


% Get average beta-weights
encoding_beta(1,1) = nanmean(glm_output.trial_type.beta_weights(1,analysis_win_idx));
encoding_beta(1,2) = nanmean(glm_output.ssd.beta_weights(1,analysis_win_idx));
encoding_beta(1,3) = nanmean(glm_output.ssd.beta_weights(2,analysis_win_idx));

end

%% Workpad: figure checks
% 
% figure;
% subplot(4,2,[1 3]); hold on
% plot(window_time,nanmean(window_sdf(behavior(behaviorIdx,:).ttx.canceled.all.all,:)),'r')
% plot(window_time,nanmean(window_sdf(behavior(behaviorIdx,:).ttx.nostop.all.all,:)),'k')
% xlim([-250 1000])
% 
% subplot(4,2,5)
% imagesc('XData',window_time,'YData',ones(1,length(window_time)),'CData',glm_output.trial_type.sig_times(1,:))
% xlim([-250 1000])
% 
% subplot(4,2,7)
% plot(window_time, glm_output.trial_type.beta_weights(1,:))
% xlim([-250 1000])
% 
% 
% subplot(4,2,[2 4]); hold on
% plot(window_time,nanmean(window_sdf(behavior(behaviorIdx,:).ttm.C.C{1},:)),'color',[1 0 0 0.33])
% plot(window_time,nanmean(window_sdf(behavior(behaviorIdx,:).ttm.C.C{2},:)),'color',[1 0 0 0.66])
% plot(window_time,nanmean(window_sdf(behavior(behaviorIdx,:).ttm.C.C{3},:)),'color',[1 0 0 0.99])
% xlim([-250 1000])
% 
% subplot(4,2,6)
% imagesc('XData',window_time,'YData',ones(1,length(window_time)),'CData',glm_output.ssd.sig_times(2,:))
% xlim([-250 1000])
% 
% subplot(4,2,8)
% plot(window_time, glm_output.ssd.beta_weights(2,:))
% xlim([-250 1000])











%% 
neuron_i = 67;

neuralFilename = mcc_map_info.session{neuron_i};
neuronLabel = mcc_map_info.unit{neuron_i};
fprintf('Extracting: %s - %s ... [%i of %i]  \n',neuralFilename,neuronLabel,neuron_i,size(mcc_map_info,1))

%... and find the corresponding behavior file index
behFilename = data_findBehFile(neuralFilename);

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

%% Curate: focus table on relevant trial types
% 


%% Setup spike data into GLM
spk_data_in = load(fullfile('D:\projects\2022-mcc-action\data\','SDF',...
    [neuralFilename '_SDF_' neuronLabel '.mat']));

event_alignment = 'stopSignal_artifical';
baseline_win = [-600:-100];

z_sdf = zscore_sdf(spk_data_in.SDF,[-600:-100],event_alignment);

[window_sdf, window_time] = movaverage_sdf(z_sdf, 100, 25);
% 2023-08-26, 23h25: I tested this and it does the job correctly. The
% z-scored profile looks exactly like the raw profile, and the movaverage
% sdf looks like a smoothed version of the raw/zscored sdf.

[n_trials, n_times] = size(window_sdf); % Get the number of windows in the time averaged data

analysis_win_idx = find(window_time >= 0 & window_time <= 600); % Find the relevant indicies for the timepoints of interest
win_fr = nanmean(window_sdf(:,analysis_win_idx),2);

reg_tbl.win_fr = win_fr; % Add the firing rate over this whole window to the GLM table

%% Run actual GLM
sig_times = []; t_betas = [];
u_sig = []; u_beta = [];

reg_tbl_canc = reg_tbl(reg_tbl.trial_type == 2,:);

% For each averaged time point
for timepoint_i = 1:n_times
    
    % Input the timepoint specific firing times
    reg_tbl_canc.firing_rate = z_sdf(reg_tbl.trial_type == 2,timepoint_i);
    
    % Fit the GLM for this
    u_t_mdl = fitlm(reg_tbl_canc,'firing_rate ~ trial_type + ssd + trial_number');
    
    % GLM output ---------------------------------
    % - Trial type
    sig_times(1,timepoint_i) = u_t_mdl.Coefficients.pValue(2) < .01; % trial type
    t_betas(1,timepoint_i) = u_t_mdl.Coefficients.tStat(2); % trial type
    
    % - Stop-signal delay
    sig_times(2,timepoint_i) = u_t_mdl.Coefficients.pValue(3) < .01; % ssd
    t_betas(2,timepoint_i) = u_t_mdl.Coefficients.tStat(3); % ssd
    
    % - Trial number
    sig_times(3,timepoint_i) = u_t_mdl.Coefficients.pValue(4) < .01; % trial_number
    t_betas(3,timepoint_i) = u_t_mdl.Coefficients.tStat(4); % trial_number

    
end


figure;
plot(window_time, t_betas(2,:))








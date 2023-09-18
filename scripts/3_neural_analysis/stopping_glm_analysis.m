
%% Neuron GLM 
% Approx 2hr run time; output saved 2023-09-03, 20h57 (stopping_glm_out)

parfor neuron_i = 1:size(mcc_map_info,1)
    neuralFilename = mcc_map_info.session{neuron_i};
    neuronLabel = mcc_map_info.unit{neuron_i};
    fprintf('Extracting: %s - %s ... [%i of %i]  \n',neuralFilename,neuronLabel,neuron_i,size(mcc_map_info,1))
    warning off
    [glm_out_mcc{neuron_i,1}, encoding_flag_mcc(neuron_i,:), encoding_beta_mcc(neuron_i,:)] =...
        stopping_glm_func(neuron_i, mcc_map_info, dataFiles_beh, dirs);
end


%% Extract stop-related spike density functions
parfor neuron_i = 1:size(mcc_map_info,1)
    
    neuralFilename = mcc_map_info.session{neuron_i};
    neuronLabel = mcc_map_info.unit{neuron_i};
    fprintf('Extracting: %s - %s ... [%i of %i]  \n',neuralFilename,neuronLabel,neuron_i,size(mcc_map_info,1))
    
    %... and find the corresponding behavior file index
    behFilename = data_findBehFile(neuralFilename);
    behaviorIdx = find(strcmp(dataFiles_beh,behFilename(1:end-4)));
    
    % Load in pre-processed spike data
    data_in = load(fullfile('D:\projects\2022-mcc-action\data\','SDF',...
        [neuralFilename '_SDF_' neuronLabel '.mat']));
    
    event_alignment = 'stopSignal_artifical';
    baseline_win = [-600:-100];
    
    z_sdf = zscore_sdf(data_in.SDF,[-600:-100],event_alignment);


    sdf_canceled_ssdx = []; sdf_nostop_ssdx = []; sdf_noncanc_ssdx = [];

    % Latency match
    for ssd_i = 1:length(behavior(behaviorIdx,:).stopSignalBeh.inh_SSD)
        trl_canceled = []; trl_canceled = behavior(behaviorIdx,:).ttm.C.C{ssd_i};
        trl_nostop = []; trl_nostop = behavior(behaviorIdx,:).ttm.C.GO{ssd_i};
        trl_noncanc = []; trl_noncanc = behavior(behaviorIdx,:).ttm.NC.NC{ssd_i};
        
        if length(trl_canceled) < 5 | length(trl_nostop) < 5 | length(trl_noncanc) < 2
            sdf_canceled_ssdx(ssd_i,:) = nan(1,length(-1000:2000));
            sdf_nostop_ssdx(ssd_i,:) = nan(1,length(-1000:2000));
            sdf_noncanc_ssdx(ssd_i,:) = nan(1,length(-1000:2000));

        else
            sdf_canceled_ssdx(ssd_i,:) = nanmean(z_sdf(trl_canceled,:));
            sdf_nostop_ssdx(ssd_i,:) = nanmean(z_sdf(trl_nostop,:));      
            sdf_noncanc_ssdx(ssd_i,:) = nanmean(z_sdf(trl_noncanc,:));      
        end
    end
    
    
    sdf_canceled{neuron_i,:} = sdf_canceled_ssdx;
    sdf_nostop{neuron_i,:} = sdf_nostop_ssdx;
    sdf_noncanc{neuron_i,:} = sdf_noncanc_ssdx;
 
    sdf_highRwd(neuron_i,:) = nanmean(z_sdf(behavior(behaviorIdx,:).ttx.canceled.all.hi,:));
    sdf_lowRwd(neuron_i,:) = nanmean(z_sdf(behavior(behaviorIdx,:).ttx.canceled.all.lo,:));
    
end

%% Export: export table for future analyses
mcc_analysis_table = table();
mcc_analysis_table.index = [1:size(mcc_map_info,1)]';
mcc_analysis_table.session = mcc_map_info.session;
mcc_analysis_table.unit = mcc_map_info.unit;

mcc_analysis_table.sdf_canceled = sdf_canceled;
mcc_analysis_table.sdf_nostop = sdf_nostop;
mcc_analysis_table.sdf_noncanc = sdf_noncanc;

mcc_analysis_table.sdf_highRwd = sdf_highRwd;
mcc_analysis_table.sdf_lowRwd = sdf_lowRwd;

mcc_analysis_table.glm_trial = encoding_flag_mcc(:,1);
mcc_analysis_table.glm_ssd_canc = encoding_flag_mcc(:,2);
mcc_analysis_table.glm_value_canc = encoding_flag_mcc(:,3);

%% Figure: produce a venn diagram of the counts (canceled)
glm_counts_canc = ...
    [... % One factor
    sum(mcc_analysis_table.glm_trial == 1 & mcc_analysis_table.glm_ssd_canc == 0 & mcc_analysis_table.glm_value_canc == 0);...
    sum(mcc_analysis_table.glm_trial == 0 & mcc_analysis_table.glm_ssd_canc == 1 & mcc_analysis_table.glm_value_canc == 0);...
    sum(mcc_analysis_table.glm_trial == 0 & mcc_analysis_table.glm_ssd_canc == 0 & mcc_analysis_table.glm_value_canc == 1);...
    % Two factors
    sum(mcc_analysis_table.glm_trial == 1 & mcc_analysis_table.glm_ssd_canc == 1 & mcc_analysis_table.glm_value_canc == 0);...
    sum(mcc_analysis_table.glm_trial == 1 & mcc_analysis_table.glm_ssd_canc == 0 & mcc_analysis_table.glm_value_canc == 1);...
    sum(mcc_analysis_table.glm_trial == 0 & mcc_analysis_table.glm_ssd_canc == 1 & mcc_analysis_table.glm_value_canc == 1);...
    % All factors
    sum(mcc_analysis_table.glm_trial == 1 & mcc_analysis_table.glm_ssd_canc == 1 & mcc_analysis_table.glm_value_canc == 1)];

p_active_neurons = ...
    (sum(mcc_analysis_table.glm_trial == 1 | mcc_analysis_table.glm_ssd_canc == 1 | mcc_analysis_table.glm_value_canc == 1))./size(mcc_analysis_table,1)

mysets = ["Action" "SSD" "Value"];
mylabels = glm_counts_canc;
venn(3,'sets',mysets,'edgeC',[0 0 0],'colors',autumn(3),'labels',mylabels,'edgeW',2);
sum(glm_counts_canc);



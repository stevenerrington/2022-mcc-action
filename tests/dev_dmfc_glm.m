
%% Neuron GLM 
% Approx 2hr run time; output saved 2023-09-03, 20h57 (stopping_glm_out)

parfor neuron_i = 1:size(dmfc_map_info,1)
    neuralFilename = dmfc_map_info.session{neuron_i};
    neuronLabel = dmfc_map_info.unit{neuron_i};
    fprintf('Extracting: %s - %s ... [%i of %i]  \n',neuralFilename,neuronLabel,neuron_i,size(dmfc_map_info,1))
    warning off
    [glm_out_dmfc{neuron_i,1}, encoding_flag_dmfc(neuron_i,:), encoding_beta_dmfc(neuron_i,:)] =...
        stopping_glm_func(neuron_i, dmfc_map_info, dataFiles_beh, dirs);
end


%% Extract stop-related spike density functions
parfor neuron_i = 1:size(dmfc_map_info,1)
    
    neuralFilename = dmfc_map_info.session{neuron_i};
    neuronLabel = dmfc_map_info.unit{neuron_i};
    fprintf('Extracting: %s - %s ... [%i of %i]  \n',neuralFilename,neuronLabel,neuron_i,size(dmfc_map_info,1))
    
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
dmfc_analysis_table = table();
dmfc_analysis_table.index = [1:size(dmfc_map_info,1)]';
dmfc_analysis_table.session = dmfc_map_info.session;
dmfc_analysis_table.unit = dmfc_map_info.unit;

dmfc_analysis_table.sdf_canceled = sdf_canceled;
dmfc_analysis_table.sdf_nostop = sdf_nostop;
dmfc_analysis_table.sdf_noncanc = sdf_noncanc;

dmfc_analysis_table.sdf_highRwd = sdf_highRwd;
dmfc_analysis_table.sdf_lowRwd = sdf_lowRwd;

dmfc_analysis_table.glm_trial = encoding_flag_dmfc(:,1);
dmfc_analysis_table.glm_ssd_canc = encoding_flag_dmfc(:,2);
dmfc_analysis_table.glm_value_canc = encoding_flag_dmfc(:,3);

%% Index GLM +ve neurons
trial_type_neurons = [];
trial_type_neurons = find(dmfc_analysis_table.glm_trial == 1);

%% Find modulation direction
glm_times = [0:10:600];
glm_index_ref = [96:156];
modulation_window = [0:600];

for neuron_i = 1:size(trial_type_neurons,1)
    neuron_j = trial_type_neurons(neuron_i);
    
    glm_window = [];
    glm_window =...
        glm_times(find(glm_out_dmfc{neuron_j}.trial_type.sig_times(1,glm_index_ref) == 1,1,'first')):...
        glm_times(find(glm_out_dmfc{neuron_j}.trial_type.sig_times(1,glm_index_ref) == 1,1,'last'));
    
    
    canc_nostop_fr_diff(neuron_i) = (nanmean(nanmean(dmfc_analysis_table.sdf_canceled{neuron_j}(:,glm_window+1000)))./...
        nanmean(nanmean(dmfc_analysis_table.sdf_nostop{neuron_j}(:,glm_window+1000))));
    
end

% Plot histogram of % diff
figure;
histogram(canc_nostop_fr_diff,0:0.05:5,'LineStyle','None'); xlim([0 2]); vline(1)
xlabel('Canceled firing % relative to no-stop'); ylabel('Count')

%% Place neurons into structure divided by type
neuron_index.trial_type.nostop = trial_type_neurons(canc_nostop_fr_diff  <= 0.9);
neuron_index.trial_type.canceled = trial_type_neurons(canc_nostop_fr_diff  >= 1.1);

dmfc_analysis_table.canc_type = ismember(1:size(dmfc_analysis_table,1),neuron_index.trial_type.canceled)';
dmfc_analysis_table.nostop_type = ismember(1:size(dmfc_analysis_table,1),neuron_index.trial_type.nostop)';

%% Output/print SDF's
for trial_type_i = {'nostop','canceled'}
    
    input_neurons = [];
    input_neurons = neuron_index.trial_type.(trial_type_i{1});
    
    n_plot_x = 4; n_plot_y = 3; n_plot_sheet = n_plot_x*n_plot_y;
    n_batches = round(size(input_neurons,1)/n_plot_sheet,-1)+1;
    
    neuron_i = 0;
    for page_i = 1:n_batches
        fig_out = figure('Renderer', 'painters', 'Position', [100 100 1200 800]);
        
        for plot_i = 1:n_plot_sheet
            neuron_i = neuron_i+1;
            try
                neuron_j = input_neurons(neuron_i);
                
                neuralFilename = dmfc_map_info.session{neuron_j};
                neuronLabel = dmfc_map_info.unit{neuron_j};
                
                %... and find the corresponding behavior file index
                behFilename = data_findBehFile(neuralFilename);
                behaviorIdx = find(strcmp(dataFiles_beh,behFilename(1:end-4)));
                
                
                subplot(n_plot_x, n_plot_y, plot_i); hold on
                plot(-1000:2000, nanmean(dmfc_analysis_table.sdf_canceled{neuron_j}),'color',colors.canceled)
                plot(-1000:2000, nanmean(dmfc_analysis_table.sdf_nostop{neuron_j}),'color',colors.nostop)
                %             plot(-1000:2000, nanmean(dmfc_analysis_table.sdf_noncanc{neuron_j}),'color',colors.noncanc)
                xlim([-200 1000]); vline(0,'k-')
                xlabel('Time from stop-signal (ms)')
                ylabel('Firing rate (spks/sec)')
                title(['Neuron: ' int2str(neuron_j) ])
                vline(behavior.stopSignalBeh(behaviorIdx).ssrt.integrationWeighted,'k--')
            catch
                continue
            end
            
        end
        
        filename = fullfile(dirs.root,'results','glm_trial_dmfc',['glm_trial_dmfc_' trial_type_i{1} '_' int2str(page_i) '_sdf.pdf']);
        set(fig_out,'PaperSize',[20 10]); %set the paper size to what you want
        print(fig_out,filename,'-dpdf') % then print it
        close(fig_out)
    end
end

%% Figure: produce a venn diagram of the counts (canceled)
dmfc_analysis_table_canceled = [];
dmfc_analysis_table_canceled = dmfc_analysis_table(neuron_index.trial_type.canceled,:);

glm_counts_canc = [];
glm_counts_canc = ...
    [... % One factor
    sum(dmfc_analysis_table_canceled.glm_trial == 1 & dmfc_analysis_table_canceled.glm_ssd_canc == 0 & dmfc_analysis_table_canceled.glm_value_canc == 0);...
    sum(dmfc_analysis_table_canceled.glm_trial == 0 & dmfc_analysis_table_canceled.glm_ssd_canc == 1 & dmfc_analysis_table_canceled.glm_value_canc == 0);...
    sum(dmfc_analysis_table_canceled.glm_trial == 0 & dmfc_analysis_table_canceled.glm_ssd_canc == 0 & dmfc_analysis_table_canceled.glm_value_canc == 1);...
    % Two factors
    sum(dmfc_analysis_table_canceled.glm_trial == 1 & dmfc_analysis_table_canceled.glm_ssd_canc == 1 & dmfc_analysis_table_canceled.glm_value_canc == 0);...
    sum(dmfc_analysis_table_canceled.glm_trial == 1 & dmfc_analysis_table_canceled.glm_ssd_canc == 0 & dmfc_analysis_table_canceled.glm_value_canc == 1);...
    sum(dmfc_analysis_table_canceled.glm_trial == 0 & dmfc_analysis_table_canceled.glm_ssd_canc == 1 & dmfc_analysis_table_canceled.glm_value_canc == 1);...
    % All factors
    sum(dmfc_analysis_table_canceled.glm_trial == 1 & dmfc_analysis_table_canceled.glm_ssd_canc == 1 & dmfc_analysis_table_canceled.glm_value_canc == 1)];


mysets = ["Action" "SSD" "Value"];
mylabels = glm_counts_canc;
venn(3,'sets',mysets,'edgeC',[0 0 0],'colors',autumn(3),'labels',mylabels,'edgeW',2);
sum(glm_counts_canc);

%% Figure: produce a venn diagram of the counts (nostop)
dmfc_analysis_table_nostop = [];
dmfc_analysis_table_nostop = dmfc_analysis_table(neuron_index.trial_type.nostop,:);

glm_counts_nostop = [];
glm_counts_nostop = ...
    [... % One factor
    sum(dmfc_analysis_table_nostop.glm_trial == 1 & dmfc_analysis_table_nostop.glm_ssd_canc == 0 & dmfc_analysis_table_nostop.glm_value_canc == 0);...
    sum(dmfc_analysis_table_nostop.glm_trial == 0 & dmfc_analysis_table_nostop.glm_ssd_canc == 1 & dmfc_analysis_table_nostop.glm_value_canc == 0);...
    sum(dmfc_analysis_table_nostop.glm_trial == 0 & dmfc_analysis_table_nostop.glm_ssd_canc == 0 & dmfc_analysis_table_nostop.glm_value_canc == 1);...
    % Two factors
    sum(dmfc_analysis_table_nostop.glm_trial == 1 & dmfc_analysis_table_nostop.glm_ssd_canc == 1 & dmfc_analysis_table_nostop.glm_value_canc == 0);...
    sum(dmfc_analysis_table_nostop.glm_trial == 1 & dmfc_analysis_table_nostop.glm_ssd_canc == 0 & dmfc_analysis_table_nostop.glm_value_canc == 1);...
    sum(dmfc_analysis_table_nostop.glm_trial == 0 & dmfc_analysis_table_nostop.glm_ssd_canc == 1 & dmfc_analysis_table_nostop.glm_value_canc == 1);...
    % All factors
    sum(dmfc_analysis_table_nostop.glm_trial == 1 & dmfc_analysis_table_nostop.glm_ssd_canc == 1 & dmfc_analysis_table_nostop.glm_value_canc == 1)];

mysets = ["Action" "SSD" "Value"];
mylabels = glm_counts_nostop;
venn(3,'sets',mysets,'edgeC',[0 0 0],'colors',summer(3),'labels',mylabels,'edgeW',2);
sum(glm_counts_nostop);



timewins.sdf = [-1000:2000];
timewins.zero = abs(timewins.sdf(1));
timewins.error_baseline = [-100:0];
timewins.error_comp = [100:300];

%% Extract: Get latency-matched SDF for non-noncanc trials
parfor neuron_i = 1:size(acc_map_info,1)
    
    neuralFilename = acc_map_info.session{neuron_i};
    neuronLabel = acc_map_info.unit{neuron_i};
    fprintf('Extracting: %s - %s ... [%i of %i]  \n',neuralFilename,neuronLabel,neuron_i,size(acc_map_info,1))
    
    %... and find the corresponding behavior file index
    behFilename = data_findBehFile(neuralFilename);
    behaviorIdx = find(strcmp(dataFiles_beh,behFilename(1:end-4)));
    
    % Load in pre-processed spike data
    data_in = load(fullfile(dirs.root,'data','SDF',...
        [neuralFilename '_SDF_' neuronLabel '.mat']));
    
    sdf_noncanc_targetx = []; sdf_nostop_targetx = [];
    sdf_noncanc_saccade = []; sdf_nostop_saccade = [];
    sdf_noncanc_tonex = []; sdf_nostop_tonex = [];
    
    % Latency match
    for ssd_i = 1:length(behavior(behaviorIdx).stopSignalBeh.inh_SSD)
        trl_noncanc = []; trl_noncanc = behavior(behaviorIdx).ttm.NC.NC{ssd_i};
        trl_nostop = []; trl_nostop = behavior(behaviorIdx).ttm.NC.GO{ssd_i};
        
        if length(trl_noncanc) < 10 | length(trl_nostop) < 10
            sdf_noncanc_saccade(ssd_i,:) = nan(1,length(-1000:2000));
            sdf_nostop_saccade(ssd_i,:) = nan(1,length(-1000:2000));
        else
            % Target:
            sdf_noncanc_targetx(ssd_i,:) = nanmean(data_in.SDF.target(trl_noncanc,:));
            sdf_nostop_targetx(ssd_i,:) = nanmean(data_in.SDF.target(trl_nostop,:));
            % SSD:
            sdf_noncanc_saccade(ssd_i,:) = nanmean(data_in.SDF.saccade(trl_noncanc,:));
            sdf_nostop_saccade(ssd_i,:) = nanmean(data_in.SDF.saccade(trl_nostop,:));
            % Tone:
            sdf_noncanc_tonex(ssd_i,:) = nanmean(data_in.SDF.tone(trl_noncanc,:));
            sdf_nostop_tonex(ssd_i,:) = nanmean(data_in.SDF.tone(trl_nostop,:));
        end
    end
    
    % Get mean SDF for:
    % - Target
    sdf_noncanc_all_target(neuron_i,:) = nanmean(sdf_noncanc_targetx);
    sdf_nostop_all_target(neuron_i,:) = nanmean(sdf_nostop_targetx);
    % - Stop-signal
    sdf_noncanc_all_saccade(neuron_i,:) = nanmean(sdf_noncanc_saccade);
    sdf_nostop_all_saccade(neuron_i,:) = nanmean(sdf_nostop_saccade);
    % - Tone
    sdf_noncanc_all_tone(neuron_i,:) = nanmean(sdf_noncanc_tonex);
    sdf_nostop_all_tone(neuron_i,:) = nanmean(sdf_nostop_tonex);

    % - Trial history
    sdf_baseline_post_nc(neuron_i,:) = nanmean(data_in.SDF.target(behavior(behaviorIdx).ttx_history.NS_after_NC,:));
    sdf_baseline_post_c(neuron_i,:) = nanmean(data_in.SDF.target(behavior(behaviorIdx).ttx_history.NS_after_C,:));
    sdf_baseline_post_ns(neuron_i,:) = nanmean(data_in.SDF.target(behavior(behaviorIdx).ttx_history.NS_after_NS,:));
    
    
    % Save SSD specific activity, just incase.
    sdf_noncanc_ssd{neuron_i} = sdf_noncanc_saccade;
    sdf_nostop_ssd{neuron_i} = sdf_nostop_saccade;
    
end

%% Analyse: Compare activity for modulation post-stopping
parfor neuron_i = 1:size(acc_map_info,1)
    
    neuralFilename = acc_map_info.session{neuron_i};
    neuronLabel = acc_map_info.unit{neuron_i};
    fprintf('Extracting: %s - %s ... [%i of %i]  \n',neuralFilename,neuronLabel,neuron_i,size(acc_map_info,1))
    
    %... and find the corresponding behavior file index
    behFilename = data_findBehFile(neuralFilename);
    behaviorIdx = find(strcmp(dataFiles_beh,behFilename(1:end-4)));
    
    % Load in pre-processed spike data
    data_in = load(fullfile(dirs.root,'data','SDF',...
        [neuralFilename '_SDF_' neuronLabel '.mat']));
    
    trl_noncanc = []; trl_nostop = [];
    trl_noncanc = behavior(behaviorIdx).ttx.noncanceled.all.all;
    trl_nostop = behavior(behaviorIdx).ttx.nostop.all.all;
    
    noncanc_sdf_mean = []; nostop_sdf_mean = [];
    noncanc_sdf_mean = nanmean(data_in.SDF.saccade(trl_noncanc,timewins.error_comp+timewins.zero),2);
    nostop_sdf_mean = nanmean(data_in.SDF.saccade(trl_nostop,timewins.error_comp+timewins.zero),2);
    
    [h,p,~,stats] = ttest2(noncanc_sdf_mean,nostop_sdf_mean);
    mean_noncanc = nanmean(noncanc_sdf_mean); mean_nostop = nanmean(nostop_sdf_mean);
    error_dir = mean_noncanc > mean_nostop;
    
    error_table(neuron_i,:) = table({neuralFilename},{neuronLabel},h,stats,error_dir,mean_noncanc,mean_nostop,...
        'VariableNames',{'neuralFilename','unit','h','stats','nc_dir','mean_nc','mean_ns'});
    
end

error_neurons_dMCC = find(error_table.h == 1 & strcmp(acc_map_info.area,'dMCC'));
error_neurons_vMCC = find(error_table.h == 1 & strcmp(acc_map_info.area,'vMCC'));

%% Analyse: Clustering approach for errors
sdfWindow = timewins.sdf;
blWindow = [-100:0];


inputSDF_target_dMCC = {sdf_noncanc_all_target(error_neurons_dMCC,:),sdf_nostop_all_target(error_neurons_dMCC,:)};
inputSDF_target_vMCC = {sdf_noncanc_all_target(error_neurons_vMCC,:),sdf_nostop_all_target(error_neurons_vMCC,:)};

inputSDF_error_dMCC = {sdf_noncanc_all_saccade(error_neurons_dMCC,:),sdf_nostop_all_saccade(error_neurons_dMCC,:)};
inputSDF_error_vMCC = {sdf_noncanc_all_saccade(error_neurons_vMCC,:),sdf_nostop_all_saccade(error_neurons_vMCC,:)};

inputSDF_tone_dMCC = {sdf_noncanc_all_tone(error_neurons_dMCC,:),sdf_nostop_all_tone(error_neurons_dMCC,:)};
inputSDF_tone_vMCC = {sdf_noncanc_all_tone(error_neurons_vMCC,:),sdf_nostop_all_tone(error_neurons_vMCC,:)};

sdfTimes = {sdfWindow, sdfWindow};
sdfEpoch = {[-200:600],[-200:600]};

colorMapping = [1,2];

%% Normalize SDFs
normResp_target_dMCC = scaleResp(inputSDF_target_dMCC,sdfTimes,'max','-bl',blWindow);
normResp_error_dMCC = scaleResp(inputSDF_error_dMCC,sdfTimes,'max','-bl',blWindow);
normResp_tone_dMCC = scaleResp(inputSDF_tone_dMCC,sdfTimes,'max','-bl',blWindow);

normResp_target_vMCC = scaleResp(inputSDF_target_vMCC,sdfTimes,'max','-bl',blWindow);
normResp_error_vMCC = scaleResp(inputSDF_error_vMCC,sdfTimes,'max','-bl',blWindow);
normResp_tone_vMCC = scaleResp(inputSDF_tone_vMCC,sdfTimes,'max','-bl',blWindow);


%% Dorsal clustering

[sortIDs_dMCC,idxDist_dMCC, raw_dMCC, respSumStruct_dMCC, rawLink_dMCC,myK_dMCC] =...
    consensusCluster(inputSDF_error_dMCC,sdfTimes,'-e',sdfEpoch,'-ei',colorMapping);

nClusters_manual = 4; clusterNeurons_dMCC = [];
for cluster_i = 1:nClusters_manual
    clusterNeurons_dMCC{cluster_i} = find(sortIDs_dMCC(:,nClusters_manual) == cluster_i );
end

% Figure: Dendrogram %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure('Renderer', 'painters', 'Position', [100 100 500 400]);

subplot(1,5,5)
[h,~,outPerm] = dendrogram(rawLink_dMCC,0,'Orientation','right');
set(gca,'YDir','Reverse');
klDendroClustChange(h,rawLink_dMCC,sortIDs_dMCC(:,nClusters_manual))
set(gca,'YTick',[]); xlabel('Similarity')

subplot(1,5,[1:4]);
for ir = 1:size(raw_dMCC,1)
    for ic = (ir+1):size(raw_dMCC,2)
        raw_dMCC(ic,ir) = raw_dMCC(ir,ic);
    end
end
imagesc(raw_dMCC(outPerm,outPerm));
colormap(gray);
xlabel('Unit Number'); set(gca,'YAxisLocation','Left');
xticks([0:100:end]); yticks([0:100:end])

% Figure: Cluster population SDF
norm_sdf_list = {normResp_target_dMCC,normResp_error_dMCC,normResp_tone_dMCC};
xlim_list = {[-250 500],[-250 500],[-500 250]};
ylim_list = {[-1 0.5],[-0.5 1],[-0.5 1],[-1 0.5]};
% Cluster main
figure('Renderer', 'painters', 'Position', [100 100 1200 900]);hold on
for epoch_i = 1:3
    clear norm_in
    norm_in = norm_sdf_list{epoch_i};
    
    for cluster_i = 1:nClusters_manual
        subplot_pos = epoch_i+(3*(cluster_i-1));
        subplot(nClusters_manual,3,subplot_pos); hold on
        plot(sdfTimes{1},nanmean(norm_in{1}(clusterNeurons_dMCC{cluster_i},:),1), 'color', [colors.noncanc]);
        plot(sdfTimes{2},nanmean(norm_in{2}(clusterNeurons_dMCC{cluster_i},:),1), 'color', [colors.nostop]);
        vline(0, 'k--'); xlim(xlim_list{epoch_i});%ylim(ylim_list{cluster_i})
        
        title(['Cluster ' int2str(cluster_i) ' - n: ' int2str(length(clusterNeurons_dMCC{cluster_i}))])
        
    end
end


%% Ventral clustering
[sortIDs_vMCC,idxDist_vMCC, raw_vMCC, respSumStruct_vMCC, rawLink_vMCC,myK_vMCC] =...
    consensusCluster(inputSDF_error_vMCC,sdfTimes,'-e',sdfEpoch,'-ei',colorMapping);

nClusters_manual = 4; clusterNeurons_vMCC = [];
for cluster_i = 1:nClusters_manual
    clusterNeurons_vMCC{cluster_i} = find(sortIDs_vMCC(:,nClusters_manual) == cluster_i );
end

% Figure: Dendrogram %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure('Renderer', 'painters', 'Position', [100 100 500 400]);

subplot(1,5,5)
[h,~,outPerm] = dendrogram(rawLink_vMCC,0,'Orientation','right');
set(gca,'YDir','Reverse');
klDendroClustChange(h,rawLink_vMCC,sortIDs_vMCC(:,nClusters_manual))
set(gca,'YTick',[]); xlabel('Similarity')

subplot(1,5,[1:4]);
for ir = 1:size(raw_vMCC,1)
    for ic = (ir+1):size(raw_vMCC,2)
        raw_vMCC(ic,ir) = raw_vMCC(ir,ic);
    end
end
imagesc(raw_vMCC(outPerm,outPerm));
colormap(gray);
xlabel('Unit Number'); set(gca,'YAxisLocation','Left');
xticks([0:100:end]); yticks([0:100:end])

% Figure: Cluster population SDF
norm_sdf_list = {normResp_target_vMCC,normResp_error_vMCC,normResp_tone_vMCC};
xlim_list = {[-250 500],[-250 500],[-500 250]};
ylim_list = {[-1 0.5],[-0.5 1],[-0.5 1],[-1 0.5]};
% Cluster main
figure('Renderer', 'painters', 'Position', [100 100 1200 900]);hold on
for epoch_i = 1:3
    clear norm_in
    norm_in = norm_sdf_list{epoch_i};
    
    for cluster_i = 1:nClusters_manual
        subplot_pos = epoch_i+(3*(cluster_i-1));
        subplot(nClusters_manual,3,subplot_pos); hold on
        plot(sdfTimes{1},nanmean(norm_in{1}(clusterNeurons_vMCC{cluster_i},:),1), 'color', [colors.noncanc]);
        plot(sdfTimes{2},nanmean(norm_in{2}(clusterNeurons_vMCC{cluster_i},:),1), 'color', [colors.nostop]);
        vline(0, 'k--'); xlim(xlim_list{epoch_i});%ylim(ylim_list{cluster_i})
        
        title(['Cluster ' int2str(cluster_i) ' - n: ' int2str(length(clusterNeurons_vMCC{cluster_i}))])
        
    end
end


%%
close all

for cluster_i = 3
    
    cluster_units = []; cluster_units = clusterNeurons_dMCC{cluster_i};
    
    cluster_map_index = error_neurons_dMCC(clusterNeurons_dMCC{cluster_i});
    
    for unit_i = 1:length(cluster_units)
        figure('Renderer', 'painters', 'Position', [100 100 700 400]);hold on
        plot(sdfTimes{1},nanmean(normResp_error_dMCC{1}(cluster_units(unit_i),:),1), 'color', [colors.noncanc]);
        plot(sdfTimes{1},nanmean(normResp_error_dMCC{2}(cluster_units(unit_i),:),1), 'color', [colors.nostop]);
        xlim([-250 1500]); vline(0,'k'); vline(600,'k'); vline(600+500,'k')
        title(['Neuron: ' int2str(cluster_map_index(unit_i)) '; Cluster ' int2str(cluster_i) ' - n: ' int2str(length(clusterNeurons_dMCC{cluster_i}))...
            ', AP: ' int2str(acc_map_info.ap(cluster_map_index(unit_i))) '; ML: ' ...
            int2str(acc_map_info.ml(cluster_map_index(unit_i))), ...
            '; Area: ' acc_map_info.area{cluster_map_index(unit_i)} ])
        
    end
end


%%


for cluster_i = 3
    
    cluster_units = []; cluster_units = clusterNeurons_vMCC{cluster_i};
    
    cluster_map_index = error_neurons_vMCC(clusterNeurons_vMCC{cluster_i});
    
    for unit_i = 1:length(cluster_units)
        figure('Renderer', 'painters', 'Position', [100 100 700 400]);hold on
        plot(sdfTimes{1},nanmean(normResp_error_vMCC{1}(cluster_units(unit_i),:),1), 'color', [colors.noncanc]);
        plot(sdfTimes{1},nanmean(normResp_error_vMCC{2}(cluster_units(unit_i),:),1), 'color', [colors.nostop]);
        xlim([-250 1500]); vline(0,'k'); vline(600,'k'); vline(600+500,'k')
        title(['Neuron: ' int2str(cluster_map_index(unit_i)) '; Cluster ' int2str(cluster_i) ' - n: ' int2str(length(clusterNeurons_dMCC{cluster_i}))...
            ', AP: ' int2str(acc_map_info.ap(cluster_map_index(unit_i))) '; ML: ' ...
            int2str(acc_map_info.ml(cluster_map_index(unit_i))), ...
            '; Area: ' acc_map_info.area{cluster_map_index(unit_i)} ])
        
    end
    
    
    
end




%%
neuron_example_list = [1309 277 1203 625];
ylims_list = {[0 10], [5 20], [0 15], [0 15]}; 

[acc_map_info.area(neuron_example_list), acc_map_info.session(neuron_example_list)]


for example_i = 1:length(neuron_example_list)
    neuron_example_index = neuron_example_list(example_i);
    
    neuralFilename = acc_map_info.session{neuron_example_index};
    neuronLabel = acc_map_info.unit{neuron_example_index};
    fprintf('Extracting: %s - %s ... [%i of %i]  \n',neuralFilename,neuronLabel,neuron_example_index,size(acc_map_info,1))
    
    %... and find the corresponding behavior file index
    behFilename = data_findBehFile(neuralFilename);
    behaviorIdx = find(strcmp(dataFiles_beh,behFilename(1:end-4)));
    
    % Load in pre-processed spike data
    data_in_sdf = load(fullfile(dirs.root,'data','SDF',...
        [neuralFilename '_SDF_' neuronLabel '.mat']));
    
    data_in_spks = load(fullfile(dirs.root,'data','Spikes',...
        [neuralFilename '_Spikes_' neuronLabel '.mat']));
    
    
    sdf_noncanc_saccade = []; sdf_nostop_saccade = [];
    
    
    % Latency match
    for ssd_i = 1:length(behavior(behaviorIdx).stopSignalBeh.inh_SSD)
        trl_noncanc = []; trl_noncanc = behavior(behaviorIdx).ttm.NC.NC{ssd_i};
        trl_nostop = []; trl_nostop = behavior(behaviorIdx).ttm.NC.GO{ssd_i};
        
        if length(trl_noncanc) < 10 | length(trl_nostop) < 10
            sdf_noncanc_saccade(ssd_i,:) = nan(1,length(-1000:2000));
            sdf_nostop_saccade(ssd_i,:) = nan(1,length(-1000:2000));
        else
            % SSD:
            sdf_noncanc_saccade(ssd_i,:) = nanmean(data_in_sdf.SDF.saccade(trl_noncanc,:));
            sdf_nostop_saccade(ssd_i,:) = nanmean(data_in_sdf.SDF.saccade(trl_nostop,:));
            
        end
    end
    
    clear noncanc_spiketimes nostop_spiketimes labels_spiketimes noncanc_sdf nostop_sdf
    
    noncanc_spiketimes = data_in_spks.Spikes.saccade(behavior(behaviorIdx).ttx.noncanceled.all.all);
    nostop_spiketimes = data_in_spks.Spikes.saccade(behavior(behaviorIdx).ttx.nostop.all.all);
    labels_spiketimes = [repmat({'1_noncanc'},length(noncanc_spiketimes),1);repmat({'2_nostop'},length(nostop_spiketimes),1)];
    
    noncanc_sdf = num2cell(data_in_sdf.SDF.saccade(behavior(behaviorIdx).ttx.noncanceled.all.all,:),2);
    nostop_sdf = num2cell(data_in_sdf.SDF.saccade(behavior(behaviorIdx).ttx.nostop.all.all,:),2);
    % noncanc_sdf = num2cell(sdf_noncanc_saccade,2);
    % nostop_sdf = num2cell(sdf_nostop_saccade,2);
    labels_sdf = [repmat({'1_noncanc'},length(noncanc_sdf),1);repmat({'2_nostop'},length(nostop_sdf),1)];
    
    clear example_sdf
    example_sdf(1,1)=gramm('x',[noncanc_spiketimes;nostop_spiketimes],'color',labels_spiketimes);
    example_sdf(1,1).geom_raster();
    example_sdf(1,1).axe_property('XLim',[-100 600]);
    
    % Produce the SDF figure
    example_sdf(2,1)=gramm('x',timewins.sdf,'y',[noncanc_sdf;nostop_sdf],'color',labels_sdf);
    example_sdf(2,1).stat_summary();
    example_sdf(2,1).axe_property('XLim',[-100 600],'YLim',ylims_list{example_i});
    example_sdf(2,1).set_names('x','Time from Target (ms)','y','FR (spk/sec)');
    example_sdf_out = figure('Renderer', 'painters', 'Position', [100 100 500 700]);
    example_sdf.draw();
    
    
end
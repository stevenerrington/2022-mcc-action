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



error_neurons = find(error_table.h == 1);


%% Analyse: Clustering approach for errors
sdfWindow = timewins.sdf;
blWindow = [-100:0];
inputNeurons = error_neurons;
inputSDF_target = {sdf_noncanc_all_target(inputNeurons,:),sdf_nostop_all_target(inputNeurons,:)};
inputSDF_error = {sdf_noncanc_all_saccade(inputNeurons,:),sdf_nostop_all_saccade(inputNeurons,:)};
inputSDF_tone = {sdf_noncanc_all_tone(inputNeurons,:),sdf_nostop_all_tone(inputNeurons,:)};

sdfTimes = {sdfWindow, sdfWindow};
sdfEpoch = {[-200:600],[-200:600]};

colorMapping = [1,2];

[sortIDs,idxDist, raw, respSumStruct, rawLink,myK] =...
    consensusCluster(inputSDF_error,sdfTimes,'-e',sdfEpoch,'-ei',colorMapping);
normResp_target = scaleResp(inputSDF_target,sdfTimes,'max','-bl',blWindow);
normResp_error = scaleResp(inputSDF_error,sdfTimes,'max','-bl',blWindow);
normResp_tone = scaleResp(inputSDF_tone,sdfTimes,'max','-bl',blWindow);

normResp_trlhistory = scaleResp({sdf_baseline_post_nc,sdf_baseline_post_c,sdf_baseline_post_ns},{sdfWindow, sdfWindow, sdfWindow},'max','-bl',blWindow);

%% 
nClusters_manual = 4; clusterNeurons = [];
for cluster_i = 1:nClusters_manual
    clusterNeurons{cluster_i} = find(sortIDs(:,nClusters_manual) == cluster_i );
end

% Figure: Dendrogram %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure('Renderer', 'painters', 'Position', [100 100 500 400]);

subplot(1,5,5)
[h,~,outPerm] = dendrogram(rawLink,0,'Orientation','right');
set(gca,'YDir','Reverse');
klDendroClustChange(h,rawLink,sortIDs(:,nClusters_manual))
set(gca,'YTick',[]); xlabel('Similarity')

subplot(1,5,[1:4]);
for ir = 1:size(raw,1)
    for ic = (ir+1):size(raw,2)
        raw(ic,ir) = raw(ir,ic);
    end
end
imagesc(raw(outPerm,outPerm));
colormap(gray);
xlabel('Unit Number'); set(gca,'YAxisLocation','Left');
xticks([0:100:end]); yticks([0:100:end])

% Figure: Cluster population SDF
norm_sdf_list = {normResp_target,normResp_error,normResp_tone};
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
        plot(sdfTimes{1},nanmean(norm_in{1}(clusterNeurons{cluster_i},:),1), 'color', [colors.noncanc]);
        plot(sdfTimes{2},nanmean(norm_in{2}(clusterNeurons{cluster_i},:),1), 'color', [colors.nostop]);
        vline(0, 'k--'); xlim(xlim_list{epoch_i});%ylim(ylim_list{cluster_i})
        
        title(['Cluster ' int2str(cluster_i) ' - n: ' int2str(length(clusterNeurons{cluster_i}))])
        
    end
end


%%
figure('Renderer', 'painters', 'Position', [100 100 600 600]);hold on

for cluster_i = 1:nClusters_manual
    subplot(nClusters_manual,1,cluster_i); hold on
    plot(sdfTimes{1},nanmean(normResp_trlhistory{1}(clusterNeurons{cluster_i},:),1), 'color', [colors.noncanc]);
    plot(sdfTimes{2},nanmean(normResp_trlhistory{2}(clusterNeurons{cluster_i},:),1), 'color', [colors.canceled]);
    plot(sdfTimes{2},nanmean(normResp_trlhistory{3}(clusterNeurons{cluster_i},:),1), 'color', [colors.nostop]);
    vline(0, 'k--'); xlim([-600 200]);
    
    title(['Cluster ' int2str(cluster_i) ' - n: ' int2str(length(clusterNeurons{cluster_i}))])
    
end




%%
close all

for cluster_i = 6
    
    cluster_units = []; cluster_units = clusterNeurons{cluster_i};
    
    for unit_i = 1:length(cluster_units)
    figure('Renderer', 'painters', 'Position', [100 100 700 400]);hold on
    plot(sdfTimes{1},nanmean(normResp_error{1}(cluster_units(unit_i),:),1), 'color', [colors.noncanc]);
    plot(sdfTimes{1},nanmean(normResp_error{2}(cluster_units(unit_i),:),1), 'color', [colors.nostop]);
    xlim([-250 1500]); vline(0,'k'); vline(600,'k'); vline(600+500,'k')
    title(['Neuron: ' int2str(cluster_units(unit_i)) '; Cluster ' int2str(cluster_i) ' - n: ' int2str(length(clusterNeurons{cluster_i}))...
        ', AP: ' int2str(acc_map_info.ap(cluster_neurons(unit_i))) '; ML: ' ...
        int2str(acc_map_info.ml(cluster_neurons(unit_i))), ...
        '; Area: ' acc_map_info.area{cluster_neurons(unit_i)} ])    
    
    end
    

    
end














%% Analyse: Spike width
spk_width_cutoff = 250;
figure('Renderer', 'painters', 'Position', [100 100 200 800]);hold on

for cluster_i = 1:nClusters_manual
    subplot(nClusters_manual,1,cluster_i)
    donut([sum(abs(acc_map_info.spk_width(inputNeurons(clusterNeurons{cluster_i}))) >= spk_width_cutoff),...
        sum(abs(acc_map_info.spk_width(inputNeurons(clusterNeurons{cluster_i}))) < spk_width_cutoff)]);
end

%% Analyse: dorsal and ventral MCC

ml_axis_range = [-8:1:8]; ap_axis_range = [20:40]; depth_range = [-3000:50:3000];
ml_axis_zero = abs(min(ml_axis_range)); ap_axis_zero = abs(min(ap_axis_range)); depth_zero = abs(min(depth_range));

figure('Renderer', 'painters', 'Position', [100 100 1400 300]);hold on
for cluster_i = 1:nClusters_manual
    cluster_neurons = []; cluster_neurons = inputNeurons(clusterNeurons{cluster_i});
    
    clear sample_coords
    sample_array = zeros(length(ml_axis_range),length(ap_axis_range), length(depth_range));
    
    for neuron_i = 1:length(cluster_neurons)
        ap_sample_ref = find(ap_axis_range == acc_map_info.ap(cluster_neurons(neuron_i)));
        ml_sample_ref = find(ml_axis_range == acc_map_info.ml(cluster_neurons(neuron_i)));
        depth_sample_ref = find(depth_range == acc_map_info.depth(cluster_neurons(neuron_i)));
        
        sample_array(ap_sample_ref,ml_sample_ref,depth_sample_ref) =...
            sample_array(ap_sample_ref,ml_sample_ref,depth_sample_ref) + 1;
        
        sample_coords(neuron_i,:) = ...
            [acc_map_info.ap(cluster_neurons(neuron_i)),...
            acc_map_info.ml(cluster_neurons(neuron_i)),...
            acc_map_info.depth(cluster_neurons(neuron_i))];
    end
    
    subplot(1,4,cluster_i); hold on
    scatter3(sample_coords(:,1),sample_coords(:,2),sample_coords(:,3),'filled')
    view(71,2.5); ylim([-8 8]); xlim([25 35]); zlim([-3000 3000])
    set(gca,'ZDir','reverse')
    xl = xlim; yl = ylim;
    patch([[1 1]*xl(1) [1 1]*xl(2)], [yl fliplr(yl)], [1 1 1 1]*0, 'k', 'FaceAlpha',0.1);
    
    grid on
end

%% Analyse: activity as a function of SSD
count = 0;
for cluster_i = 1:nClusters_manual
    neuron_list = [];
    neuron_list = inputNeurons(clusterNeurons{cluster_i});
    
    for neuron_list_i = 1:length(neuron_list)
        neuron_i = neuron_list(neuron_list_i);
        count = count + 1;
        behFilename = data_findBehFile(acc_map_info.session{neuron_list(neuron_list_i)});
        behaviorIdx = find(strcmp(dataFiles_beh,behFilename(1:end-4)));
        
        mid_ssd_idx = behavior(behaviorIdx).stopSignalBeh.midSSDidx;
        
        conflict_epoch = [0:100]+timewins.zero;
        
        early_ssd_activity = mean(sdf_noncanc_ssd{neuron_i}(mid_ssd_idx-1,conflict_epoch));
        mid_ssd_activity = mean(sdf_noncanc_ssd{neuron_i}(mid_ssd_idx,conflict_epoch));
        late_ssd_activity = mean(sdf_noncanc_ssd{neuron_i}(mid_ssd_idx+1,conflict_epoch));
        
        all_ssd_activity_norm = [early_ssd_activity;mid_ssd_activity;late_ssd_activity] ./...
            max([early_ssd_activity;mid_ssd_activity;late_ssd_activity]);
        
        
        early_ssd_ssd = behavior(behaviorIdx).stopSignalBeh.inh_SSD(mid_ssd_idx-1);
        mid_ssd_ssd = behavior(behaviorIdx).stopSignalBeh.inh_SSD(mid_ssd_idx);
        late_ssd_ssd = behavior(behaviorIdx).stopSignalBeh.inh_SSD(mid_ssd_idx+1);
        
        early_ssd_pnc = behavior(behaviorIdx).stopSignalBeh.inh_pnc(mid_ssd_idx-1);
        mid_ssd_pnc = behavior(behaviorIdx).stopSignalBeh.inh_pnc(mid_ssd_idx);
        late_ssd_pnc = behavior(behaviorIdx).stopSignalBeh.inh_pnc(mid_ssd_idx+1);
        
        
        functional_stopping_activity(count,:) =...
            table(neuron_i,cluster_i,...
            all_ssd_activity_norm(1), all_ssd_activity_norm(2), all_ssd_activity_norm(3),...
            early_ssd_ssd,mid_ssd_ssd,late_ssd_ssd,early_ssd_pnc,mid_ssd_pnc,late_ssd_pnc,...
            'VariableNames',{'neuron','cluster','activity_early','activity_mid','activity_late',...
            'ssd_early','ssd_mid','ssd_late','pnc_early','pnc_mid','pnc_late'});
        
    end
end

%% Figure: Conflict Activity
cluster_i = 1;

plotData = []; plotDataY = []; labelData = [];

plotData=...
    [functional_stopping_activity.activity_early;...
    functional_stopping_activity.activity_mid;...
    functional_stopping_activity.activity_late];

plotDataY=...
    [functional_stopping_activity.pnc_early;...
    functional_stopping_activity.pnc_mid;...
    functional_stopping_activity.pnc_late];

labelData=...
    [repmat({'1_early'},length(functional_stopping_activity.activity_early),1);...
    repmat({'2_mid'},length(functional_stopping_activity.activity_mid),1);...
    repmat({'3_late'},length(functional_stopping_activity.activity_late),1)];

clusterData = repmat(functional_stopping_activity.cluster,3,1);

testFigure_out = figure('Renderer', 'painters', 'Position', [100 100 1000 200]);
clear testFigure
testFigure(1,1)= gramm('x',labelData,'y',plotData,'color',labelData);
testFigure(1,1).stat_summary('geom',{'bar','line','black_errorbar'});
testFigure.axe_property('YLim',[0.8 1]);
testFigure.facet_grid([],clusterData);
testFigure.draw();

testFigure_out2 = figure('Renderer', 'painters', 'Position', [100 100 1000 200]);
clear testFigure2
testFigure2(1,1)= gramm('x',plotDataY,'y',plotData);
testFigure2(1,1).stat_glm();
testFigure2(1,1).geom_point();
testFigure2.facet_grid([],clusterData);
testFigure2.draw();


%% PCA Analysis: Error
for neuron_i = 1:size(acc_map_info,1)
    
    neuralFilename = acc_map_info.session{neuron_i};
    neuronLabel = acc_map_info.unit{neuron_i};    
    %... and find the corresponding behavior file index
    behFilename = data_findBehFile(neuralFilename);
    behaviorIdx = find(strcmp(dataFiles_beh,behFilename(1:end-4)));
    mid_ssd_idx = behavior(behaviorIdx).stopSignalBeh.midSSDidx;
 
    nc_ssd1(neuron_i,:) = sdf_noncanc_ssd{neuron_i}(mid_ssd_idx-1,:);
    nc_ssd2(neuron_i,:) = sdf_noncanc_ssd{neuron_i}(mid_ssd_idx,:);
    nc_ssd3(neuron_i,:) = sdf_noncanc_ssd{neuron_i}(mid_ssd_idx+1,:);
    
    ns_ssd1(neuron_i,:) = sdf_nostop_ssd{neuron_i}(mid_ssd_idx-1,:);
    ns_ssd2(neuron_i,:) = sdf_nostop_ssd{neuron_i}(mid_ssd_idx,:);
    ns_ssd3(neuron_i,:) = sdf_nostop_ssd{neuron_i}(mid_ssd_idx+1,:);   
end

% Set data parameters & windows
ephysWin = [-1000 2000]; winZero = abs(ephysWin(1));
plotWin = [-250:500]; 
analyseTime = [-100:200];
getColors

% Set PCA parameters
samplingRate = 1/1000; % inherent to the data. Do not change
numPCs = 8; % pick a number that will capture most of the variance
timeStep = 10; % smaller values will yield high resolution at the expense of computation time, default will sample at 20ms
withinConditionsOnly = false; % if true, will only analyze tanlging for times within the same condition

clear PCA_mainInput PCA_mainOutput
% Setup data for PCA analysis
PCA_mainInput(1).A = sdf_noncanc_all_saccade';
PCA_mainInput(1).times = plotWin';
PCA_mainInput(1).analyzeTimes = analyseTime';
PCA_mainInput(1).condition = 'SSD2 - Error';
PCA_mainInput(1).color = [colors.noncanc 1];

PCA_mainInput(2).A = sdf_nostop_all_saccade';
PCA_mainInput(2).times = plotWin';
PCA_mainInput(2).analyzeTimes = analyseTime';
PCA_mainInput(2).condition = 'SSD2 - Correct';
PCA_mainInput(2).color = [colors.nostop 1];


% Run PCA analysis
[~, PCA_mainOutput] = tangleAnalysis(PCA_mainInput, samplingRate,'numPCs',numPCs,'softenNorm',5 ,'timeStep',timeStep,'withinConditionsOnly',withinConditionsOnly); % soft normalize neural data
pca_figure(PCA_mainInput,PCA_mainOutput)

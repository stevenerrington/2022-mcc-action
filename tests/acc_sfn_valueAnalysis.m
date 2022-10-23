timewins.sdf = [-1000:2000];
timewins.zero = abs(timewins.sdf(1));
timewins.reward_baseline = [-100:0];
timewins.reward_comp = [0:250];

%% Extract: Get latency-matched SDF for non-canceled trials
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
 
    % Get mean SDF for:
    % - Target
    sdf_hiRwd_target(neuron_i,:) = nanmean(data_in.SDF.target(behavior(behaviorIdx).ttx.nostop.all.hi,:));
    sdf_loRwd_target(neuron_i,:) = nanmean(data_in.SDF.target(behavior(behaviorIdx).ttx.nostop.all.lo,:));
    % - Stop-signal
    sdf_hiRwd_saccade(neuron_i,:) = nanmean(data_in.SDF.saccade(behavior(behaviorIdx).ttx.nostop.all.hi,:));
    sdf_loRwd_saccade(neuron_i,:) = nanmean(data_in.SDF.saccade(behavior(behaviorIdx).ttx.nostop.all.lo,:));
    % - SSRT
    sdf_hiRwd_tone(neuron_i,:) =  nanmean(data_in.SDF.tone(behavior(behaviorIdx).ttx.nostop.all.hi,:));
    sdf_loRwd_tone(neuron_i,:) =  nanmean(data_in.SDF.tone(behavior(behaviorIdx).ttx.nostop.all.lo,:));
    % - Tone
    sdf_hiRwd_reward(neuron_i,:) =  nanmean(data_in.SDF.reward(behavior(behaviorIdx).ttx.nostop.all.hi,:));
    sdf_loRwd_reward(neuron_i,:) =  nanmean(data_in.SDF.reward(behavior(behaviorIdx).ttx.nostop.all.lo,:));
    
    
    [h,p,~,stats] = ttest2(...
        nanmean(data_in.SDF.reward(behavior(behaviorIdx).ttx.nostop.all.lo,timewins.reward_comp+timewins.zero),2),...
        nanmean(data_in.SDF.reward(behavior(behaviorIdx).ttx.nostop.all.hi,timewins.reward_comp+timewins.zero),2));
    mean_hi = nanmean(nanmean(data_in.SDF.reward(behavior(behaviorIdx).ttx.nostop.all.hi,timewins.reward_comp+timewins.zero),2));
    mean_lo = nanmean(nanmean(data_in.SDF.reward(behavior(behaviorIdx).ttx.nostop.all.lo,timewins.reward_comp+timewins.zero),2));
    rwd_dir = mean_hi > mean_lo;
    
    rwd_table(neuron_i,:) = table({neuralFilename},{neuronLabel},h,stats,rwd_dir,mean_hi,mean_lo,...
        'VariableNames',{'neuralFilename','unit','h','stats','rwd_dir','mean_hi','mean_lo'});
    
    
end

rwd_neurons = find(rwd_table.h == 1);


%% Analyse: Clustering approach for errors
sdfWindow = timewins.sdf;
blWindow = timewins.reward_baseline;
inputNeurons = rwd_neurons;
inputSDF_target = {sdf_loRwd_target(inputNeurons,:),sdf_hiRwd_target(inputNeurons,:)};
inputSDF_saccade = {sdf_loRwd_saccade(inputNeurons,:),sdf_hiRwd_saccade(inputNeurons,:)};
inputSDF_tone = {sdf_loRwd_tone(inputNeurons,:),sdf_hiRwd_tone(inputNeurons,:)};
inputSDF_reward = {sdf_loRwd_reward(inputNeurons,:),sdf_hiRwd_reward(inputNeurons,:)};

sdfTimes = {sdfWindow, sdfWindow};
sdfEpoch = {[-200:600],[-200:600]};

colorMapping = [1,2];

[sortIDs,idxDist, raw, respSumStruct, rawLink,myK] =...
    consensusCluster(inputSDF_reward,sdfTimes,'-e',sdfEpoch,'-ei',colorMapping);
normResp_target = scaleResp(inputSDF_target,sdfTimes,'max','-bl',blWindow);
normResp_saccade = scaleResp(inputSDF_saccade,sdfTimes,'max','-bl',blWindow);
normResp_tone = scaleResp(inputSDF_tone,sdfTimes,'max','-bl',blWindow);
normResp_reward = scaleResp(inputSDF_reward,sdfTimes,'max','-bl',blWindow);

nClusters_manual = 6; clusterNeurons = [];
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
norm_sdf_list = {normResp_target,normResp_saccade,normResp_tone,normResp_reward};
xlim_list = {[-250 500],[-500 500],[-250 500],[-500 500]};
ylim_list = {[-0.5 1.5],[-0.5 1.5],[-1 0.5],[-1 0.5]};
% Cluster main
figure('Renderer', 'painters', 'Position', [100 100 1200 900]);hold on
for epoch_i = 1:4
    clear norm_in
    norm_in = norm_sdf_list{epoch_i};
    
    for cluster_i = 1:nClusters_manual
        subplot_pos = epoch_i+(4*(cluster_i-1));
        subplot(nClusters_manual,4,subplot_pos); hold on
        plot(sdfTimes{1},nanmean(norm_in{1}(clusterNeurons{cluster_i},:),1), 'color', [0,188,227]/255);
        plot(sdfTimes{2},nanmean(norm_in{2}(clusterNeurons{cluster_i},:),1), 'color', [255, 0, 255]/255);
        vline(0, 'k--'); xlim(xlim_list{epoch_i});%ylim(ylim_list{cluster_i})
        
        title(['Cluster ' int2str(cluster_i) ' - n: ' int2str(length(clusterNeurons{cluster_i}))])
        
    end
end

%% Analyse: Spike width
spk_width_cutoff = 250;
figure('Renderer', 'painters', 'Position', [100 100 200 800]);hold on

for cluster_i = 1:nClusters_manual
    subplot(nClusters_manual,1,cluster_i)
    pie([sum(abs(acc_map_info.spk_width(inputNeurons(clusterNeurons{cluster_i}))) >= spk_width_cutoff),...
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
    
    subplot(1,nClusters_manual,cluster_i); hold on
    scatter3(sample_coords(:,1),sample_coords(:,2),sample_coords(:,3),'filled')
    view(71,2.5); ylim([-8 8]); xlim([25 35]); zlim([-3000 3000])
    set(gca,'ZDir','reverse')
    xl = xlim; yl = ylim;
    patch([[1 1]*xl(1) [1 1]*xl(2)], [yl fliplr(yl)], [1 1 1 1]*0, 'k', 'FaceAlpha',0.1);
    
    grid on
end


%% PCA Analysis: Reward - target
% Set data parameters & windows
ephysWin = [-1000 2000]; winZero = abs(ephysWin(1));
plotWin = [-500:1000]; 
analyseTime = [-500:1000];
getColors

% Set PCA parameters
samplingRate = 1/1000; % inherent to the data. Do not change
numPCs = 8; % pick a number that will capture most of the variance
timeStep = 10; % smaller values will yield high resolution at the expense of computation time, default will sample at 20ms
withinConditionsOnly = false; % if true, will only analyze tanlging for times within the same condition

clear PCA_mainInput PCA_mainOutput
% Setup data for PCA analysis
PCA_mainInput(1).A = sdf_loRwd_target';
PCA_mainInput(1).times = plotWin';
PCA_mainInput(1).analyzeTimes = analyseTime';
PCA_mainInput(1).condition = 'Low Reward';
PCA_mainInput(1).color = [0,188,227]/255;

% Setup data for PCA analysis
PCA_mainInput(2).A = sdf_hiRwd_target';
PCA_mainInput(2).times = plotWin';
PCA_mainInput(2).analyzeTimes = analyseTime';
PCA_mainInput(2).condition = 'High Reward';
PCA_mainInput(2).color = [255, 0, 255]/255;


% Run PCA analysis
[~, PCA_mainOutput] = tangleAnalysis(PCA_mainInput, samplingRate,'numPCs',numPCs,'softenNorm',5 ,'timeStep',timeStep,'withinConditionsOnly',withinConditionsOnly); % soft normalize neural data
pca_figure(PCA_mainInput,PCA_mainOutput)

%% PCA Analysis: Reward - tone
% Set data parameters & windows
ephysWin = [-1000 2000]; winZero = abs(ephysWin(1));
plotWin = [-500:1000]; 
analyseTime = [-500:1000];
getColors

% Set PCA parameters
samplingRate = 1/1000; % inherent to the data. Do not change
numPCs = 8; % pick a number that will capture most of the variance
timeStep = 10; % smaller values will yield high resolution at the expense of computation time, default will sample at 20ms
withinConditionsOnly = false; % if true, will only analyze tanlging for times within the same condition

clear PCA_mainInput PCA_mainOutput
% Setup data for PCA analysis
PCA_mainInput(1).A = sdf_loRwd_tone';
PCA_mainInput(1).times = plotWin';
PCA_mainInput(1).analyzeTimes = analyseTime';
PCA_mainInput(1).condition = 'Low Reward';
PCA_mainInput(1).color = [0,188,227]/255;

% Setup data for PCA analysis
PCA_mainInput(2).A = sdf_hiRwd_tone';
PCA_mainInput(2).times = plotWin';
PCA_mainInput(2).analyzeTimes = analyseTime';
PCA_mainInput(2).condition = 'High Reward';
PCA_mainInput(2).color = [255, 0, 255]/255;


% Run PCA analysis
[~, PCA_mainOutput] = tangleAnalysis(PCA_mainInput, samplingRate,'numPCs',numPCs,'softenNorm',5 ,'timeStep',timeStep,'withinConditionsOnly',withinConditionsOnly); % soft normalize neural data
pca_figure(PCA_mainInput,PCA_mainOutput)

%% 
%% Extract SDF

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
    
    sdf_canceled_ssdx = []; sdf_nostop_ssdx = [];
    
    % Latency match
    for ssd_i = 1:length(behavior(behaviorIdx,:).stopSignalBeh.inh_SSD)
        trl_canceled = []; trl_canceled = behavior(behaviorIdx,:).ttm.C.C{ssd_i};
        trl_nostop = []; trl_nostop = behavior(behaviorIdx,:).ttm.C.GO{ssd_i};
        
        if length(trl_canceled) < 10 | length(trl_nostop) < 10
            sdf_canceled_ssdx(ssd_i,:) = nan(1,length(-1000:2000));
            sdf_nostop_ssdx(ssd_i,:) = nan(1,length(-1000:2000));
        else           
            sdf_canceled_ssdx(ssd_i,:) = nanmean(data_in.SDF.stopSignal_artifical(trl_canceled,:));
            sdf_nostop_ssdx(ssd_i,:) = nanmean(data_in.SDF.stopSignal_artifical(trl_nostop,:));
        end
    end
    
    norm_mean = nanmean(nanmean([sdf_canceled_ssdx, sdf_nostop_ssdx]));
    norm_std = nanstd(nanmean([sdf_canceled_ssdx, sdf_nostop_ssdx]));
    
    sdf_canceled_all_stopsignal(neuron_i,:) = (nanmean(sdf_canceled_ssdx)-norm_mean)./norm_std;
    sdf_nostop_all_stopsignal(neuron_i,:) = (nanmean(sdf_nostop_ssdx)-norm_mean)./norm_std;
    
end


%% Analyse: Clustering approach for errors

inputSDF = {sdf_canceled_all_stopsignal,sdf_nostop_all_stopsignal};

sdfTimes = {[-1000:2000],[-1000:2000]};
sdfEpoch = {[0:500],[0:500]};

colorMapping = [1,1];

[sortIDs,idxDist, raw, respSumStruct, rawLink,myK] =...
    consensusCluster(inputSDF,sdfTimes,'-e',sdfEpoch,'-ei',colorMapping,'-er',sdfEpoch,'-c',0.5);

nClusters_manual = myK; clusterNeurons = [];
for cluster_i = 1:nClusters_manual
    clusterNeurons{cluster_i} = find(sortIDs(:,nClusters_manual) == cluster_i );
end


% Generate a quick sdf of each cluster
for cluster_i = 1:nClusters_manual
    figure('Renderer', 'painters', 'Position', [100 100 500 300]);hold on

    a = subplot(1,1,1); hold on
    plot(sdfTimes{1},nanmean(sdf_canceled_all_stopsignal(clusterNeurons{cluster_i},:),1), 'color', [colors.canceled]);
    plot(sdfTimes{1},nanmean(sdf_nostop_all_stopsignal(clusterNeurons{cluster_i},:),1), 'color', [colors.nostop]);
    vline([0], 'k--'); xlim([-500 250])


    xlim([-200 1500])
    
    title(['Cluster ' int2str(cluster_i) ' - n: ' int2str(length(clusterNeurons{cluster_i}))])
    
    pause
    close all
end


% Generate a quick sdf of each cluster
for cluster_i = [5 11 16 19 23]
    figure('Renderer', 'painters', 'Position', [100 100 500 300]);hold on

    a = subplot(1,1,1); hold on
    plot(sdfTimes{1},sdf_canceled_all_stopsignal(clusterNeurons{cluster_i},:), 'color', [colors.canceled 0.1]);
    plot(sdfTimes{1},sdf_nostop_all_stopsignal(clusterNeurons{cluster_i},:), 'color', [colors.nostop 0.1]);
    
    plot(sdfTimes{1},nanmean(sdf_canceled_all_stopsignal(clusterNeurons{cluster_i},:)), 'color', [colors.canceled],'LineWidth',2);
    plot(sdfTimes{1},nanmean(sdf_nostop_all_stopsignal(clusterNeurons{cluster_i},:)), 'color', [colors.nostop],'LineWidth',2);
    
    vline([0], 'k--'); xlim([-500 250])


    xlim([-200 1500])
    
    title(['Cluster ' int2str(cluster_i) ' - n: ' int2str(length(clusterNeurons{cluster_i}))])
    
    pause
    close all
end



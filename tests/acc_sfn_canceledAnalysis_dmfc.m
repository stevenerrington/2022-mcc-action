

dajo_datamap_curated_DMFC = data_sessionCurate...
    (dajo_datamap,...
    'area', {'DMFC'}, 'monkey', {'dar','jou'}, 'signal', {'SPK'}, 'spacing', [50, 100, 150]);



% Get session information and indices
n_sessions = size(dajo_datamap_curated_DMFC,1);
dmfc_map_info = [];

dataFiles_beh = unique(dajo_datamap_curated_DMFC.sessionBeh);
dataFiles_neural = unique(dajo_datamap_curated_DMFC.dataFilename);

behavior = acc_stopping_extractBeh(dirs,dataFiles_beh);



for session_i = 1:n_sessions
    neuralFilename = dajo_datamap_curated_DMFC.dataFilename{session_i};
    behFilename = data_findBehFile(neuralFilename);
    beh_index = util_find_beh_index(behavior,behFilename);
    logInfo(session_i,:) = util_getLogInfo(neuralFilename);
    
    fprintf(['Extracting data for ' neuralFilename ': session %i of %i     \n'],...
        session_i, n_sessions);
        
    spk_data = load(fullfile(dirs.data,[neuralFilename '-spk.mat']));
    
    n_neurons= size(dajo_datamap_curated_DMFC.spkInfo(session_i,1).unitInfo,1);
    map_info = table();
    
    for neuron_i = 1:n_neurons
        site = dajo_datamap_curated_DMFC.spkInfo(session_i,1).unitInfo.site(neuron_i);
        session = dajo_datamap_curated_DMFC.dataFilename(session_i);
        monkey = dajo_datamap_curated_DMFC.monkey(session_i);
        unit = dajo_datamap_curated_DMFC.spkInfo(session_i,1).unitInfo.unitDSP(neuron_i);
        ap = logInfo.ap_stereo(session_i);
        ml = logInfo.ml_stereo(session_i);
        mua = dajo_datamap_curated_DMFC.spkInfo(session_i,1).unitInfo.flag_mua(neuron_i);
        noise = dajo_datamap_curated_DMFC.spkInfo(session_i,1).unitInfo.flag_noise(neuron_i);

        spk_width = util_getSpkWidth(spk_data,dajo_datamap_curated_DMFC.spkInfo(session_i,1).unitInfo.unitWAV{neuron_i});
        
        map_info(neuron_i,:) = table(session,monkey,unit,site,ap,ml,mua,noise,spk_width);
    end
    
    dmfc_map_info = [dmfc_map_info; map_info];
end

% Clear up unknown areas and remove noise clusters
dmfc_map_info(dmfc_map_info.noise == 1,:) = [];


%%
spk_width_cutoff = 250;

% 
% %%
% dmfc_session_list = unique(dmfc_map_info.session);
% 
% parfor (session_i = 1:length(dmfc_session_list),4)
% 
%     % We get the (neural) filename of the record of interest
%     neuralFilename = dmfc_session_list{session_i};
%     fprintf('Extracting: %s ... [%i of %i]  \n',neuralFilename,length(dmfc_session_list))
%     
%     %... and find the corresponding behavior file index
%     behFilename = data_findBehFile(neuralFilename);
%     behaviorIdx = find(strcmp(dataFiles_beh,behFilename(1:end-4)));
%     
%     import_data = struct(); import_data = load_spkFile(dirs,neuralFilename);
%     
%     % Convolve spike times to get continous trace
%     spk_data_sdf = [];spk_data_spikes = [];
%     [spk_data_sdf, spk_data_spikes] =...
%         spk_alignTrials(behavior(behaviorIdx).trialEventTimes(:,[3,5,6,7,9,10]),...
%         import_data.time, [-1000 2000]);
% 
%     % Then split this data into individual channels
%     neuron_labels = {};
%     neuron_labels = dmfc_map_info.unit(strcmp(dmfc_map_info.session,neuralFilename));
%     
%     for neuron_i = 1:length(neuron_labels)
%         SDF = []; Spikes = [];
%         SDF = spk_data_sdf.(neuron_labels{neuron_i});
%         Spikes = spk_data_spikes.(neuron_labels{neuron_i});
%         
%         save_filename_sdf = ...
%             fullfile(dirs.root,'data','SDF',...
%             [neuralFilename '_SDF_' neuron_labels{neuron_i} '.mat']);
%         save_filename_spk = ...
%             fullfile(dirs.root,'data','Spikes',...
%             [neuralFilename '_Spikes_' neuron_labels{neuron_i} '.mat']);
%         
%         util_parsaveSDF(save_filename_sdf, SDF)
%         util_parsaveSPK(save_filename_spk, Spikes)
%     end
%     
% end
% 

%% Extract: Get latency-matched SDF for non-canceled trials
timewins.sdf = [-1000:2000];
timewins.zero = abs(timewins.sdf(1));
timewins.ssrt_baseline = [-100:0];
timewins.ssrt_comp = [-50:250];

parfor neuron_i = 1:size(dmfc_map_info,1)
    
    neuralFilename = dmfc_map_info.session{neuron_i};
    neuronLabel = dmfc_map_info.unit{neuron_i};
    fprintf('Extracting: %s - %s ... [%i of %i]  \n',neuralFilename,neuronLabel,neuron_i,size(dmfc_map_info,1))
    
    %... and find the corresponding behavior file index
    behFilename = data_findBehFile(neuralFilename);
    behaviorIdx = find(strcmp(dataFiles_beh,behFilename(1:end-4)));
    
    % Load in pre-processed spike data
    data_in = load(fullfile(dirs.root,'data','SDF',...
        [neuralFilename '_SDF_' neuronLabel '.mat']));
    
    sdf_canceled_targetx = []; sdf_nostop_targetx = [];
    sdf_canceled_ssdx = []; sdf_nostop_ssdx = [];
    sdf_canceled_ssrtx = []; sdf_nostop_ssrtx = [];
    sdf_canceled_tonex = []; sdf_nostop_tonex = [];
    
    % Latency match
    for ssd_i = 1:length(behavior(behaviorIdx).stopSignalBeh.inh_SSD)
        trl_canceled = []; trl_canceled = behavior(behaviorIdx).ttm.C.C{ssd_i};
        trl_nostop = []; trl_nostop = behavior(behaviorIdx).ttm.C.GO{ssd_i};
        
        if length(trl_canceled) < 10 | length(trl_nostop) < 10
            sdf_canceled_ssdx(ssd_i,:) = nan(1,length(-1000:2000));
            sdf_nostop_ssdx(ssd_i,:) = nan(1,length(-1000:2000));
        else
            % Target:
            sdf_canceled_targetx(ssd_i,:) = nanmean(data_in.SDF.target(trl_canceled,:));
            sdf_nostop_targetx(ssd_i,:) = nanmean(data_in.SDF.target(trl_nostop,:));
            % SSD:
            sdf_canceled_ssdx(ssd_i,:) = nanmean(data_in.SDF.stopSignal_artifical(trl_canceled,:));
            sdf_nostop_ssdx(ssd_i,:) = nanmean(data_in.SDF.stopSignal_artifical(trl_nostop,:));
            % SSRT:
            sdf_canceled_ssrtx(ssd_i,:) = nanmean(data_in.SDF.ssrt(trl_canceled,:));
            sdf_nostop_ssrtx(ssd_i,:) = nanmean(data_in.SDF.ssrt(trl_nostop,:));
            % Tone:
            sdf_canceled_tonex(ssd_i,:) = nanmean(data_in.SDF.tone(trl_canceled,:));
            sdf_nostop_tonex(ssd_i,:) = nanmean(data_in.SDF.tone(trl_nostop,:));
        end
    end
    
    trl_noncanceled = [];
    trl_noncanceled = behavior(behaviorIdx).ttx.noncanceled.all.all;
    
    % Get mean SDF for:
    % - Target
    sdf_canceled_all_target(neuron_i,:) = nanmean(sdf_canceled_targetx);
    sdf_nostop_all_target(neuron_i,:) = nanmean(sdf_nostop_targetx);
    % - Stop-signal
    sdf_canceled_all_stopsignal(neuron_i,:) = nanmean(sdf_canceled_ssdx);
    sdf_nostop_all_stopsignal(neuron_i,:) = nanmean(sdf_nostop_ssdx);
    sdf_noncanc_all_stopsignal(neuron_i,:) = nanmean(data_in.SDF.stopSignal_artifical(trl_noncanceled,:));
    % - SSRT
    sdf_canceled_all_ssrt(neuron_i,:) = nanmean(sdf_canceled_ssrtx);
    sdf_nostop_all_ssrt(neuron_i,:) = nanmean(sdf_nostop_ssrtx);
    % - Tone
    sdf_canceled_all_tone(neuron_i,:) = nanmean(sdf_canceled_tonex);
    sdf_nostop_all_tone(neuron_i,:) = nanmean(sdf_nostop_tonex);
    sdf_noncanc_all_tone(neuron_i,:) = nanmean(data_in.SDF.tone(trl_noncanceled,:));
    
     % Error ---- IN PROGRESS
    sdf_nostop_saccade(neuron_i,:) = nanmean(data_in.SDF.saccade(behavior(behaviorIdx).ttx.nostop.all.all,:));
    sdf_noncanc_saccade(neuron_i,:) = nanmean(data_in.SDF.saccade(behavior(behaviorIdx).ttx.noncanceled.all.all,:));
        
    
    % Value ---- IN PROGRESS
    sdf_canceled_hi_stopsignal(neuron_i,:) = nanmean(data_in.SDF.stopSignal_artifical(behavior(behaviorIdx).ttx.canceled.all.hi,:));
    sdf_canceled_lo_stopsignal(neuron_i,:) = nanmean(data_in.SDF.stopSignal_artifical(behavior(behaviorIdx).ttx.canceled.all.lo,:));
    sdf_nostop_hi_stopsignal(neuron_i,:) = nanmean(data_in.SDF.stopSignal_artifical(behavior(behaviorIdx).ttx.nostop.all.hi,:));
    sdf_nostop_lo_stopsignal(neuron_i,:) = nanmean(data_in.SDF.stopSignal_artifical(behavior(behaviorIdx).ttx.nostop.all.lo,:));
    
    % Trial history ---- IN PROGRESS
    sdf_post_canc_target(neuron_i,:) = nanmean(data_in.SDF.target(behavior(behaviorIdx).ttx_history.NS_after_C,:));
    sdf_post_nostop_target(neuron_i,:) = nanmean(data_in.SDF.target(behavior(behaviorIdx).ttx_history.NS_after_NS,:));
    sdf_post_noncanc_target(neuron_i,:) = nanmean(data_in.SDF.target(behavior(behaviorIdx).ttx_history.NS_after_NC,:));

    % Save SSD specific activity, just incase.
    sdf_canceled_ssd{neuron_i} = sdf_canceled_ssdx;
    sdf_nostop_ssd{neuron_i} = sdf_nostop_ssdx;
    
end

%% Analyse: Compare activity for modulation post-stopping
parfor neuron_i = 1:size(dmfc_map_info,1)
    
    neuralFilename = dmfc_map_info.session{neuron_i};
    neuronLabel = dmfc_map_info.unit{neuron_i};
    fprintf('Extracting: %s - %s ... [%i of %i]  \n',neuralFilename,neuronLabel,neuron_i,size(dmfc_map_info,1))
    
    %... and find the corresponding behavior file index
    behFilename = data_findBehFile(neuralFilename);
    behaviorIdx = find(strcmp(dataFiles_beh,behFilename(1:end-4)));
    
    % Load in pre-processed spike data
    data_in = load(fullfile(dirs.root,'data','SDF',...
        [neuralFilename '_SDF_' neuronLabel '.mat']));
    
    trl_canceled = []; trl_nostop = [];
    trl_canceled = behavior(behaviorIdx).ttx.canceled.all.all;
    trl_nostop = behavior(behaviorIdx).ttx.nostop.all.all;
    
    canc_sdf_mean = []; nostop_sdf_mean = [];
    canc_sdf_mean = nanmean(data_in.SDF.ssrt(trl_canceled,timewins.ssrt_comp+timewins.zero),2);
    nostop_sdf_mean = nanmean(data_in.SDF.ssrt(trl_nostop,timewins.ssrt_comp+timewins.zero),2);
    
    [h,p,~,stats] = ttest2(canc_sdf_mean,nostop_sdf_mean);
    mean_canceled = nanmean(canc_sdf_mean); mean_nostop = nanmean(nostop_sdf_mean);
    ssrt_dir = mean_canceled > mean_nostop;
    
    ssrt_table(neuron_i,:) = table({neuralFilename},{neuronLabel},h,stats,ssrt_dir,mean_canceled,mean_nostop,...
        'VariableNames',{'neuralFilename','unit','h','stats','c_dir','mean_c','mean_ns'});
    
end



ssrt_neurons = find(ssrt_table.h == 1);



%% Analyse: Clustering approach for stopping
sdfWindow = timewins.sdf;
blWindow = [-100:0];
inputNeurons = ssrt_neurons;
inputSDF_target = {sdf_canceled_all_target(inputNeurons,:),sdf_nostop_all_target(inputNeurons,:)};
inputSDF_ssd = {sdf_canceled_all_stopsignal(inputNeurons,:),sdf_nostop_all_stopsignal(inputNeurons,:),sdf_noncanc_all_stopsignal(inputNeurons,:)};
inputSDF_ssrt = {sdf_canceled_all_ssrt(inputNeurons,:),sdf_nostop_all_ssrt(inputNeurons,:)};
inputSDF_tone = {sdf_canceled_all_tone(inputNeurons,:),sdf_nostop_all_tone(inputNeurons,:),sdf_noncanc_all_tone(inputNeurons,:)};

sdfTimes = {sdfWindow, sdfWindow};
sdfEpoch = {[200:600],[200:600]};

colorMapping = [1,2];

[sortIDs,idxDist, raw, respSumStruct, rawLink,myK] =...
    consensusCluster(inputSDF_ssrt,sdfTimes,'-e',sdfEpoch,'-ei',colorMapping);
normResp_target = scaleResp(inputSDF_target,sdfTimes,'max','-bl',blWindow);
normResp_stopSignal = scaleResp(inputSDF_ssd,[sdfTimes, timewins.sdf],'max','-bl',blWindow);
normResp_ssrt = scaleResp(inputSDF_ssrt,sdfTimes,'max','-bl',blWindow);
normResp_tone = scaleResp(inputSDF_tone,[sdfTimes, timewins.sdf],'max','-bl',blWindow);

normResp_error = scaleResp(...
    {sdf_nostop_saccade(inputNeurons,:),sdf_noncanc_saccade(inputNeurons,:)},...
    repmat({timewins.sdf},2,1),'max','-bl',blWindow);

normResp_value = scaleResp(...
    {sdf_canceled_lo_stopsignal(inputNeurons,:),sdf_canceled_hi_stopsignal(inputNeurons,:),...
    sdf_nostop_lo_stopsignal(inputNeurons,:),sdf_nostop_hi_stopsignal(inputNeurons,:)},...
    repmat({timewins.sdf},4,1),'max','-bl',blWindow);

normResp_trialhistory = scaleResp(...
    {sdf_post_canc_target(inputNeurons,:),sdf_post_nostop_target(inputNeurons,:),...
    sdf_post_noncanc_target(inputNeurons,:)},repmat({timewins.sdf},3,1),'max','-bl',blWindow);

nClusters_manual = 10; clusterNeurons = [];
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
norm_sdf_list = {normResp_target,normResp_stopSignal,normResp_ssrt,normResp_tone};
xlim_list = {[-250 500],[-250 500],[-250 500],[-500 1000]};
ylim_list = {[-0.5 1.5],[-0.5 1.5],[-1 0.5],[-1 0.5]};
% Cluster main
figure('Renderer', 'painters', 'Position', [100 100 1200 900]);hold on
for epoch_i = 1:4
    clear norm_in
    norm_in = norm_sdf_list{epoch_i};
    
    for cluster_i = 1:nClusters_manual
        subplot_pos = epoch_i+(4*(cluster_i-1));
        subplot(nClusters_manual,4,subplot_pos); hold on
        plot(sdfTimes{1},nanmean(norm_in{1}(clusterNeurons{cluster_i},:),1), 'color', [colors.canceled]);
        plot(sdfTimes{2},nanmean(norm_in{2}(clusterNeurons{cluster_i},:),1), 'color', [colors.nostop]);
        if epoch_i == 2 || epoch_i ==4 
        plot(sdfTimes{2},nanmean(norm_in{3}(clusterNeurons{cluster_i},:),1), 'color', [colors.noncanc]);
        end
        vline(0, 'k--'); xlim(xlim_list{epoch_i});%ylim(ylim_list{cluster_i})
        title(['Cluster ' int2str(cluster_i) ' - n: ' int2str(length(clusterNeurons{cluster_i}))])
        
    end
end



%% Analyse: activity as a function of SSD
count = 0;
for cluster_i = 1:nClusters_manual
    neuron_list = [];
    neuron_list = inputNeurons(clusterNeurons{cluster_i});
    
    for neuron_list_i = 1:length(neuron_list)
        neuron_i = neuron_list(neuron_list_i);
        count = count + 1;
        behFilename = data_findBehFile(dmfc_map_info.session{neuron_list(neuron_list_i)});
        behaviorIdx = find(strcmp(dataFiles_beh,behFilename(1:end-4)));
        ssrt = round(behavior(behaviorIdx).stopSignalBeh.ssrt.integrationWeighted);
        
        mid_ssd_idx = behavior(behaviorIdx).stopSignalBeh.midSSDidx;
        
        conflict_epoch = []; conflict_epoch = [-50:50]+ssrt+timewins.zero;
        
        early_ssd_activity_canc = mean(sdf_canceled_ssd{neuron_i}(mid_ssd_idx-1,conflict_epoch));
        mid_ssd_activity_canc = mean(sdf_canceled_ssd{neuron_i}(mid_ssd_idx,conflict_epoch));
        late_ssd_activity_canc = mean(sdf_canceled_ssd{neuron_i}(mid_ssd_idx+1,conflict_epoch));
        
        early_ssd_activity_nostop = mean(sdf_nostop_ssd{neuron_i}(mid_ssd_idx-1,conflict_epoch));
        mid_ssd_activity_nostop = mean(sdf_nostop_ssd{neuron_i}(mid_ssd_idx,conflict_epoch));
        late_ssd_activity_nostop = mean(sdf_nostop_ssd{neuron_i}(mid_ssd_idx+1,conflict_epoch));
        
        all_ssd_activity_norm = [early_ssd_activity_canc;mid_ssd_activity_canc;late_ssd_activity_canc]...
            ./max([early_ssd_activity_canc;mid_ssd_activity_canc;late_ssd_activity_canc]);     
        
        sdf_max = max([early_ssd_activity_canc;mid_ssd_activity_canc;late_ssd_activity_canc;...
            early_ssd_activity_nostop;mid_ssd_activity_nostop;late_ssd_activity_nostop]);
        
        ssd_sdf_all_canc{cluster_i,1}(neuron_list_i,:) = sdf_canceled_ssd{neuron_i}(mid_ssd_idx-1,:)./sdf_max;
        ssd_sdf_all_canc{cluster_i,2}(neuron_list_i,:) = sdf_canceled_ssd{neuron_i}(mid_ssd_idx,:)./sdf_max;
        ssd_sdf_all_canc{cluster_i,3}(neuron_list_i,:) = sdf_canceled_ssd{neuron_i}(mid_ssd_idx+1,:)./sdf_max;

        ssd_sdf_all_nostop{cluster_i,1}(neuron_list_i,:) = sdf_nostop_ssd{neuron_i}(mid_ssd_idx-1,:)./sdf_max;
        ssd_sdf_all_nostop{cluster_i,2}(neuron_list_i,:) = sdf_nostop_ssd{neuron_i}(mid_ssd_idx,:)./sdf_max;
        ssd_sdf_all_nostop{cluster_i,3}(neuron_list_i,:) = sdf_nostop_ssd{neuron_i}(mid_ssd_idx+1,:)./sdf_max;       
        
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

testFigure_out = figure('Renderer', 'painters', 'Position', [100 100 1000 400]);
clear testFigure
testFigure(1,1)= gramm('x',labelData,'y',plotData,'color',labelData);
testFigure(1,1).stat_summary('geom',{'point','line','errorbar'});
testFigure(2,1)= gramm('x',labelData,'y',plotDataY,'color',labelData);
testFigure(2,1).stat_summary('geom',{'point','line','errorbar'});

testFigure(1,1).axe_property('YLim',[0.8 1]);
testFigure(2,1).axe_property('YLim',[0.0 1]);
testFigure(1,1).facet_grid([],clusterData);
testFigure(2,1).facet_grid([],clusterData);
testFigure.draw();

%% Figure;

figure('Renderer', 'painters', 'Position', [100 100 1000 200]);hold on

for cluster_i = 1:nClusters_manual
    subplot(1,nClusters_manual,cluster_i); hold on
    plot(-1000:2000,nanmean(ssd_sdf_all_canc{cluster_i,1}),'color',[colors.canceled 0.3])
    plot(-1000:2000,nanmean(ssd_sdf_all_canc{cluster_i,2}),'color',[colors.canceled 0.6])
    plot(-1000:2000,nanmean(ssd_sdf_all_canc{cluster_i,3}),'color',[colors.canceled 1.0])
    
    plot(-1000:2000,nanmean(ssd_sdf_all_nostop{cluster_i,1}),'color',[colors.nostop 0.3])
    plot(-1000:2000,nanmean(ssd_sdf_all_nostop{cluster_i,2}),'color',[colors.nostop 0.6])
    plot(-1000:2000,nanmean(ssd_sdf_all_nostop{cluster_i,3}),'color',[colors.nostop 1.0])
    xlim([-200 600]); vline(0,'k')
end








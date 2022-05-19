
%% Setup workspace
clear all; clc; % Clear workspace
addpath('C:\Users\Steven\Desktop\Projects\2022-dajo-toolbox'); % DaJo-toolbox
addpath('C:\Users\Steven\Desktop\Projects\2022-acc-stopping\_toolbox\clustering_schall') % Clustering toolbox
dirs = data_setDir();

%% Curate data
dajo_datamap = load_datamap(dirs);

dajo_datamap_curated = data_sessionCurate...
    (dajo_datamap,...
    'area', {'ACC'}, 'monkey', {'dar','jou'}, 'signal', {'SPK'}, 'spacing', [50, 100, 150]);

dataFiles_beh = unique(dajo_datamap_curated.sessionBeh);
dataFiles_neural = unique(dajo_datamap_curated.dataFilename);

%% Extract behavioral data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear behavior

% Looping through each of the individual data files
parfor dataFileIdx = 1:length(dataFiles_beh)
    % We first report loop status:
    fprintf('Extracting: %s ... [%i of %i]  \n',dataFiles_beh{dataFileIdx},dataFileIdx,length(dataFiles_beh))
    
    % We then get the (behavior) filename of the record of interest
    behFilename = [dataFiles_beh{dataFileIdx} '-beh'];
    % and load it into the workspace
    import_data = struct(); import_data = load_behFile(dirs,behFilename);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Once we have the behavior in matlab, we can then look to extract
    % relevant task/session behavior.
    
    sessionName = {behFilename}; % Get the session name
    [ttx, ~, trialEventTimes] =... % Index of trials and event timings
        beh_getTrials(import_data.events.stateFlags_,import_data.events.Infos_);
    [stopSignalBeh, ~] = beh_getStoppingInfo... % Stopping behavior
        (import_data.events.stateFlags_,import_data.events.Infos_,ttx);
    trialEventTimes.ssrt = trialEventTimes.stopSignal_artifical + stopSignalBeh.ssrt.integrationWeighted;
    [ttm] = beh_getMatchedTrials(stopSignalBeh,ttx, trialEventTimes); % Trial matching indices
        
    % After extracting the individual behavioral variable, we then collapse
    % it into one structure for the given session.
    behavior(dataFileIdx) = struct('sessionName',sessionName,'ttx',ttx,'trialEventTimes',trialEventTimes,...
        'stopSignalBeh',stopSignalBeh,'ttm',ttm);
    
end


%% Extract neural data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
parfor dataFileIdx = 1:length(dataFiles_neural)
    % We first report loop status:
    fprintf('Extracting: %s ... [%i of %i]  \n',dataFiles_neural{dataFileIdx},dataFileIdx,length(dataFiles_neural))
    
    % We then get the (neural) filename of the record of interest
    neuralFilename = dataFiles_neural{dataFileIdx};
   
    %... and find the corresponding behavior file index
    behFilename = data_findBehFile(neuralFilename);
    behaviorIdx = find(strcmp(dataFiles_beh,behFilename(1:end-4)));
    
    % Load in data-file
    import_data = struct(); import_data = load_spkFile(dirs,neuralFilename);
    
    % Convolve spike times to get continous trace
    spk_data = [];
    spk_data = spk_alignTrials(behavior(behaviorIdx).trialEventTimes(:,[3,9,10]),...
        import_data.time, [-1000 2000]);

    
    ssd_in = [behavior(behaviorIdx).stopSignalBeh.midSSDidx - 1,...
        behavior(behaviorIdx).stopSignalBeh.midSSDidx,...
        behavior(behaviorIdx).stopSignalBeh.midSSDidx + 1];
    
    signal_average{dataFileIdx} = neural_ttmSignalAverage(spk_data,...
        behavior(behaviorIdx).ttm.C,...
         'units',fieldnames(spk_data),...
         'events',{'target','stopSignal_artifical','ssrt'},...
         'ssd',ssd_in,...
         'weighting',{'on'})
    
end

signal_collapse = neural_collapseSignalSession(signal_average,...
    'events',{'target','stopSignal_artifical','ssrt'},...
    'conditions',{'C','GO'},...
    'conditions_map',[1 2]);

%% Clustering

sdfWindow = [-200:600];
inputSDF = {signal_collapse.stopSignal_artifical.C(:,sdfWindow+1000),...
    signal_collapse.stopSignal_artifical.GO(:,sdfWindow+1000)};

sdfTimes = {sdfWindow, sdfWindow};
sdfEpoch = {[0:600],[0:600]};

colorMapping = [1,2];

[sortIDs,idxDist, raw, respSumStruct, rawLink,myK] =...
    consensusCluster(inputSDF,sdfTimes,'-e',sdfEpoch,'-ei',colorMapping);
normResp = scaleResp(inputSDF,sdfTimes,'max');

nClusters_manual = 5; clusterNeurons = [];
for i = 1:nClusters_manual
    clusterNeurons{i} = find(sortIDs(:,nClusters_manual) == i );
end

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



for i = 1:nClusters_manual
    figure('Renderer', 'painters', 'Position', [100 100 500 400]);hold on
    
    plot(sdfTimes{1},nanmean(normResp{1}(clusterNeurons{i},:),1), 'color', 'r');
    plot(sdfTimes{2},nanmean(normResp{2}(clusterNeurons{i},:),1), 'color', 'b');
    vline(0, 'k--'); xlim([-200 600])
    
    title(['Cluster ' int2str(i) ' - n: ' int2str(length(clusterNeurons{i}))])
    
end

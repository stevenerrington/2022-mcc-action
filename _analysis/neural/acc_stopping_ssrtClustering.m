
%% Clustering
sdfWindow = [-200:500];
inputSDF = {signal_collapse.ssrt.C(:,sdfWindow+1000),...
    signal_collapse.ssrt.GO(:,sdfWindow+1000)};

sdfTimes = {sdfWindow, sdfWindow};
sdfEpoch = {[0:500],[0:500]};

colorMapping = [1,2];

[sortIDs,idxDist, raw, respSumStruct, rawLink,myK] =...
    consensusCluster(inputSDF,sdfTimes,'-e',sdfEpoch,'-ei',colorMapping);
normResp = scaleResp(inputSDF,sdfTimes,'max','-bl',[0:500]);

nClusters_manual = 6; clusterNeurons = [];
for i = 1:nClusters_manual
    clusterNeurons{i} = find(sortIDs(:,nClusters_manual) == i );
end

%% Figure: Dendrogram
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


%% Cluster main
figure('Renderer', 'painters', 'Position', [100 100 1200 200]);hold on

for i = 1:nClusters_manual
    subplot(1,nClusters_manual,i); hold on
    plot(sdfTimes{1},nanmean(normResp{1}(clusterNeurons{i},:),1), 'color', [colors.canceled]);
    plot(sdfTimes{2},nanmean(normResp{2}(clusterNeurons{i},:),1), 'color', [colors.nostop]);
    vline(0, 'k--'); xlim([-200 600])
    
    title(['Cluster ' int2str(i) ' - n: ' int2str(length(clusterNeurons{i}))])
    
end

%% Cluster SSD
count = 0;

for session_i = 1:size(dajo_datamap_curated,1)
    for unit_i = 1:size(signal_average_spk{session_i}.individual.ssrt,1)
        count = count+1;
        sdf_ssd_ns_1(count,:) = signal_average_spk{session_i}.individual.ssrt.signal_ssd{unit_i,1}{1, 1}(1,:);
        sdf_ssd_ns_2(count,:) = signal_average_spk{session_i}.individual.ssrt.signal_ssd{unit_i,1}{1, 1}(2,:);
        sdf_ssd_ns_3(count,:) = signal_average_spk{session_i}.individual.ssrt.signal_ssd{unit_i,1}{1, 1}(3,:);
        
        sdf_ssd_c_1(count,:) = signal_average_spk{session_i}.individual.ssrt.signal_ssd{unit_i,1}{1, 2}(1,:);
        sdf_ssd_c_2(count,:) = signal_average_spk{session_i}.individual.ssrt.signal_ssd{unit_i,1}{1, 2}(2,:);
        sdf_ssd_c_3(count,:) = signal_average_spk{session_i}.individual.ssrt.signal_ssd{unit_i,1}{1, 2}(3,:);
        
    end
end

clear sdfWindow inputSDF_ssd
sdfWindow = [-200:500];
inputSDF_ssd = {sdf_ssd_c_1(:,sdfWindow+1000),...
    sdf_ssd_c_2(:,sdfWindow+1000),...
    sdf_ssd_c_3(:,sdfWindow+1000),...
    sdf_ssd_ns_1(:,sdfWindow+1000),...
    sdf_ssd_ns_2(:,sdfWindow+1000),...
    sdf_ssd_ns_3(:,sdfWindow+1000)};
sdfTimes_SSD = {sdfWindow,sdfWindow,sdfWindow,sdfWindow,sdfWindow,sdfWindow};
normResp_SSD = scaleResp(inputSDF_ssd,sdfTimes_SSD,'max','-bl',[0:500]);

figure('Renderer', 'painters', 'Position', [100 100 1200 200]);hold on
for i = 1:nClusters_manual
    subplot(1,nClusters_manual,i); hold on
    
    plot(sdfTimes{1},nanmean(normResp_SSD{1}(clusterNeurons{i},:),1), 'color', [colors.canceled 0.3]);
    plot(sdfTimes{2},nanmean(normResp_SSD{2}(clusterNeurons{i},:),1), 'color', [colors.canceled 0.6]);
    plot(sdfTimes{2},nanmean(normResp_SSD{3}(clusterNeurons{i},:),1), 'color', [colors.canceled 0.9]);
    
    plot(sdfTimes{2},nanmean(normResp_SSD{4}(clusterNeurons{i},:),1), 'color', [colors.nostop 0.3]);
    plot(sdfTimes{2},nanmean(normResp_SSD{5}(clusterNeurons{i},:),1), 'color', [colors.nostop 0.25]);
    plot(sdfTimes{2},nanmean(normResp_SSD{6}(clusterNeurons{i},:),1), 'color', [colors.nostop 0.25]);
    vline(0, 'k--'); xlim([-200 600])
    
    title(['Cluster ' int2str(i) ' - n: ' int2str(length(clusterNeurons{i}))])
end

%% Conflict slope
% Get pNC array for each session
for session_i = 1:length(behavior)

    ssd_idx = [behavior(session_i).stopSignalBeh.midSSDidx-1,...
        behavior(session_i).stopSignalBeh.midSSDidx,...
        behavior(session_i).stopSignalBeh.midSSDidx+1];
    
    beh_pnc_session(session_i,:) = ...
        behavior(session_i).stopSignalBeh.inh_pnc(ssd_idx);
        
end

pnc_beh_unit = [];
for ii = 1:size(dajo_datamap_curated,1)
    nUnits = size(dajo_datamap_curated.spkInfo(ii).unitInfo,1);
    beh_i = find(strcmp(dajo_datamap_curated.sessionBeh(ii),dataFiles_beh));
    
    pnc_beh_unit = [pnc_beh_unit; repmat(beh_pnc_session(beh_i,:),nUnits,1)];
end

conflictWindow = [0:250]+abs(sdfWindow(1));

for cond_i = 1:6
    norm_fr_all(:,cond_i)=nanmean(normResp_SSD{cond_i}(:,conflictWindow),2);
end

norm_fr_c = norm_fr_all(:,1:3);
norm_fr_ns = norm_fr_all(:,3:end);

for unit_i = 1:length(norm_fr_c)
    [pnc_fr_r(unit_i,1), pnc_fr_r(unit_i,2)] = corr(norm_fr_c(unit_i,:)',pnc_beh_unit(unit_i,:)');
end


for i = 1:nClusters_manual
    figure('Renderer', 'painters', 'Position', [100 100 500 400]);hold on
    
    histogram(pnc_fr_r(clusterNeurons{i},1),-1:0.1:1)
    
    title(['Cluster ' int2str(i) ' - n: ' int2str(length(clusterNeurons{i}))])
end




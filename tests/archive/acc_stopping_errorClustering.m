
%% Clustering
sdfWindow = [-999:2000];
blWindow = [0:100];
event = 'saccade';
inputSDF = {signal_collapse_noncanc.(event).NC(:,sdfWindow+1000),...
    signal_collapse_noncanc.(event).GO(:,sdfWindow+1000)};

sdfTimes = {sdfWindow, sdfWindow};
sdfEpoch = {[0:500],[0:500]};

colorMapping = [1,2];

[sortIDs,idxDist, raw, respSumStruct, rawLink,myK] =...
    consensusCluster(inputSDF,sdfTimes,'-e',sdfEpoch,'-ei',colorMapping);
normResp = scaleResp(inputSDF,sdfTimes,'max','-bl',blWindow);

nClusters_manual = 5; clusterNeurons = [];
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
    plot(sdfTimes{1},nanmean(normResp{1}(clusterNeurons{i},:),1), 'color', [colors.noncanc]);
    plot(sdfTimes{2},nanmean(normResp{2}(clusterNeurons{i},:),1), 'color', [colors.nostop]);
    vline(0, 'k--'); xlim([-200 600])
    
    title(['Cluster ' int2str(i) ' - n: ' int2str(length(clusterNeurons{i}))])
    
end

%% 
cluster_i = 2;

clusterNeurons{cluster_i}












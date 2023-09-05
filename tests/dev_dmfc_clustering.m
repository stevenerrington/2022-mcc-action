
%% Setup & parameterize clustering
input_neurons = [];
input_neurons = neuron_index.trial_type.canceled;

canc_sdf_in = []; nostop_sdf_in = [];

for neuron_i = 1:length(input_neurons)
    neuron_j = input_neurons(neuron_i);
    canc_sdf_in(neuron_i,:) = nanmean(dmfc_analysis_table.sdf_canceled{neuron_j});
    nostop_sdf_in(neuron_i,:) = nanmean(dmfc_analysis_table.sdf_nostop{neuron_j});
    diff_sdf_in(neuron_i,:) = canc_sdf_in(neuron_i,:)-nostop_sdf_in(neuron_i,:);
end

inputSDF = [];
inputSDF = {canc_sdf_in,nostop_sdf_in};

sdfTimes = {[-1000:2000],[-1000:2000]};
sdfEpoch = {[0:600],[0:600]};

colorMapping = [1,1];

%% Clustering algorithm
[sortIDs,idxDist, raw, respSumStruct, rawLink,idK] =...
    consensusCluster(inputSDF,sdfTimes,'-e',sdfEpoch,'-ei',colorMapping,'-er',sdfEpoch,'-c',0.5);

% Refine number of clusters (iterative with below dendrogram)
myK = 3; % idK = based on consensusCluster raw output

%% Plot dendrogram
figure('Renderer', 'painters', 'Position', [100 100 500 400]);

subplot(1,5,5)
[h,~,outPerm] = dendrogram(rawLink,0,'Orientation','right');
set(gca,'YDir','Reverse');
klDendroClustChange(h,rawLink,sortIDs(:,myK))
set(gca,'YTick',[],'YLim',[1 length(rawLink)]); xlabel('Similarity')

subplot(1,5,[1:4]);
for ir = 1:size(raw,1)
    for ic = (ir+1):size(raw,2)
        raw(ic,ir) = raw(ir,ic);
    end
end
imagesc(raw(outPerm,outPerm));
colormap(flipud(viridis));
xlabel('Unit Number'); set(gca,'YAxisLocation','Left');
set(gca,'CLim',[-1 1])


%% Identify/assign neurons to clusters
nClusters_manual = myK; clusterNeurons = [];
for cluster_i = 1:nClusters_manual
    clusterNeurons{cluster_i} = find(sortIDs(:,nClusters_manual) == cluster_i );
end

%% Plot cluster populations
% Generate a quick population sdf of each cluster
for cluster_i = 1:nClusters_manual
    figure('Renderer', 'painters', 'Position', [100 100 500 300]);hold on

    example_neuron_subfig = subplot(1,1,1); hold on
    plot(sdfTimes{1},nanmean(canc_sdf_in(clusterNeurons{cluster_i},:),1), 'color', [colors.canceled]);
    plot(sdfTimes{1},nanmean(nostop_sdf_in(clusterNeurons{cluster_i},:),1), 'color', [colors.nostop]);
    vline([0], 'k--'); xlim([-250 750])
    
    title(['Cluster ' int2str(cluster_i) ' - n: ' int2str(length(clusterNeurons{cluster_i}))])

end

% Print each neuron within a cluster
for cluster_i = 1:nClusters_manual
    
    input_neurons_cluster = [];
    input_neurons_cluster = clusterNeurons{cluster_i};
    
    n_plot_x = 4; n_plot_y = 3; n_plot_sheet = n_plot_x*n_plot_y;
    n_batches = ceil(size(input_neurons_cluster,1)/n_plot_sheet);
    
    neuron_i = 0;
    for page_i = 1:n_batches
        fig_out = figure('Renderer', 'painters', 'Position', [100 100 1200 800]);
        
        for plot_i = 1:n_plot_sheet
            neuron_i = neuron_i+1;
            try
                neuron_j = input_neurons_cluster(neuron_i);
                
                neuralFilename = dmfc_map_info.session{neuron_j};
                neuronLabel = dmfc_map_info.unit{neuron_j};
                
                %... and find the corresponding behavior file index
                behFilename = data_findBehFile(neuralFilename);
                behaviorIdx = find(strcmp(dataFiles_beh,behFilename(1:end-4)));
                
                
                subplot(n_plot_x, n_plot_y, plot_i); hold on
                plot(-1000:2000, canc_sdf_in(neuron_j,:),'color',colors.canceled)
                plot(-1000:2000, nostop_sdf_in(neuron_j,:),'color',colors.nostop)
                xlim([-250 750]); vline(0,'k-')
                xlabel('Time from stop-signal (ms)')
                ylabel('Firing rate (spks/sec)')
                title(['Cluster: ' int2str(cluster_i) ' | Neuron: ' int2str(neuron_j) ])
                vline(behavior.stopSignalBeh(behaviorIdx).ssrt.integrationWeighted,'k--')
            catch
                continue
            end
            
        end
        
        filename = fullfile(dirs.root,'results','cluster_dmfc',['cluster_dmfc_canc_' int2str(cluster_i) '_pg' int2str(page_i) '.pdf']);
        set(fig_out,'PaperSize',[20 10]); %set the paper size to what you want
        print(fig_out,filename,'-dpdf') % then print it
        close(fig_out)
    end
end

%% Heatmap: Cluster activity
% Generate a quick population sdf of each cluster
figure('Renderer', 'painters', 'Position', [100 100 750 400]);hold on

for cluster_i = 1:nClusters_manual
    sdf_in_diff = []; sdf_in_diff = diff_sdf_in(clusterNeurons{cluster_i},:);
    sdf_in_diff = sdf_in_diff./max(abs(sdf_in_diff(:,[0:600]+1000)),[],2);
    
    subplot(3,myK,[cluster_i]); hold on
    plot(sdfTimes{1},nanmean(canc_sdf_in(clusterNeurons{cluster_i},:),1), 'color', [colors.canceled]);
    plot(sdfTimes{1},nanmean(nostop_sdf_in(clusterNeurons{cluster_i},:),1), 'color', [colors.nostop]);
    vline([0], 'k--'); xlim([-250 750])
    
    subplot(3,myK,[cluster_i+myK cluster_i+(myK*2)]); hold on
    imagesc('XData',[-1000:2000],'YData',[1:size(sdf_in_diff,1)],'CData',sdf_in_diff)
    colormap('viridis')
    vline([0], 'k--'); xlim([-250 750]); ylim([1-0.5 max(cellfun(@length,clusterNeurons))+0.5]); caxis([0 1])
    title(['Cluster ' int2str(cluster_i) ' - n: ' int2str(length(clusterNeurons{cluster_i}))])
    set(gca,'YDir','Reverse')
    
end

%% Figure: Raster, example SDF, population SDF

% Example neuron idx ----------------------------
% Cluster 1: 1
% Cluster 2: 51
% Cluster 3: 73/166/246
% Cluster 4: 30
% Cluster 5: 32/105/199

example_neuron_idx = [1, 51, 166, 30, 105];
example_neuron_ylim = {[0 15],[0 15],[5 25],[0 10],[5 15]};
for cluster_i = 1:myK
    input_neuron_idx = example_neuron_idx(cluster_i);
    raw_neuron_idx = input_neurons(input_neuron_idx);
    
    main_sdf_fig_out = [];
    main_sdf_fig_out = main_sdf_figure(dirs,colors,...
        dmfc_analysis_table,raw_neuron_idx,behavior,dataFiles_beh,...
        clusterNeurons, cluster_i, canc_sdf_in, nostop_sdf_in);
    
    main_sdf_fig_out(2,1).axe_property('YLim',example_neuron_ylim{cluster_i});
    figure('Renderer', 'painters', 'Position', [100 100 300 500]);
    main_sdf_fig_out.draw
end





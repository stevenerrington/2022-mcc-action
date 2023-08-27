timewins.sdf = [-1000:2000];
timewins.zero = abs(timewins.sdf(1));
timewins.error_baseline = [-300:-100];
timewins.error_comp = [100:400];

%% Extract: Get latency-matched SDF for non-noncanc trials
parfor neuron_i = 1:size(mcc_map_info,1)
    
    neuralFilename = mcc_map_info.session{neuron_i};
    neuronLabel = mcc_map_info.unit{neuron_i};
    fprintf('Extracting: %s - %s ... [%i of %i]  \n',neuralFilename,neuronLabel,neuron_i,size(mcc_map_info,1))
    
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
    for ssd_i = 1:length(behavior.stopSignalBeh(behaviorIdx).inh_SSD)
        trl_noncanc = []; trl_noncanc = behavior.ttm(behaviorIdx).NC.NC{ssd_i};
        trl_nostop = []; trl_nostop = behavior.ttm(behaviorIdx).NC.GO{ssd_i};
        
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
    sdf_baseline_post_nc(neuron_i,:) = nanmean(data_in.SDF.target(behavior.ttx_history(behaviorIdx).NS_after_NC,:));
    sdf_baseline_post_c(neuron_i,:) = nanmean(data_in.SDF.target(behavior.ttx_history(behaviorIdx).NS_after_C,:));
    sdf_baseline_post_ns(neuron_i,:) = nanmean(data_in.SDF.target(behavior.ttx_history(behaviorIdx).NS_after_NS,:));
    
    
    % Save SSD specific activity, just incase.
    sdf_noncanc_ssd{neuron_i} = sdf_noncanc_saccade;
    sdf_nostop_ssd{neuron_i} = sdf_nostop_saccade;
    
end

%% Analyse: Get normalized activity
for neuron_i = 1:size(mcc_map_info,1)

    bl_mean = nanmean([sdf_noncanc_all_target(neuron_i,timewins.error_baseline+timewins.zero),...
        sdf_nostop_all_target(neuron_i,timewins.error_baseline+timewins.zero)]);
    bl_std = nanstd([sdf_noncanc_all_target(neuron_i,timewins.error_baseline+timewins.zero),...
        sdf_nostop_all_target(neuron_i,timewins.error_baseline+timewins.zero)]);
    

    norm_sdf.target.noncanc(neuron_i,:) = (sdf_noncanc_all_target(neuron_i,:)-bl_mean)./bl_std;
    norm_sdf.target.nostop(neuron_i,:) = (sdf_nostop_all_target(neuron_i,:)-bl_mean)./bl_std;
    norm_sdf.target.diff(neuron_i,:) = norm_sdf.target.noncanc(neuron_i,:) -  norm_sdf.target.nostop(neuron_i,:);
    
    norm_sdf.saccade.noncanc(neuron_i,:) = (sdf_noncanc_all_saccade(neuron_i,:)-bl_mean)./bl_std;
    norm_sdf.saccade.nostop(neuron_i,:) = (sdf_nostop_all_saccade(neuron_i,:)-bl_mean)./bl_std;
    norm_sdf.saccade.diff(neuron_i,:) = norm_sdf.saccade.noncanc(neuron_i,:) -  norm_sdf.saccade.nostop(neuron_i,:);

    norm_sdf.tone.noncanc(neuron_i,:) = (sdf_noncanc_all_tone(neuron_i,:)-bl_mean)./bl_std;
    norm_sdf.tone.nostop(neuron_i,:) = (sdf_nostop_all_tone(neuron_i,:)-bl_mean)./bl_std;
    norm_sdf.tone.diff(neuron_i,:) = norm_sdf.tone.noncanc(neuron_i,:) -  norm_sdf.tone.nostop(neuron_i,:);

end


%% Figure: Plot normalized activity as a heatmap
[~,sort_activity_idx] = sort(nanmean(norm_sdf.saccade.diff(:,timewins.error_comp+timewins.zero),2));

figure('Renderer', 'painters', 'Position', [100 100 300 400]);
imagesc(-1000:2000,1:size(mcc_map_info,1),norm_sdf.saccade.diff(sort_activity_idx,:))
xlim([-200 600]); set(gca,'CLim',[-3 3]); colorbar('southoutside')
set(gca,'TickDir','out','FontSize',10,'YTick',[],'YColor',[1 1 1]); box off
set(gcf,'color','w');
colormap(colors.error_gradient); vline(0,'k')


%% Analyse: Compare activity for modulation post-stopping
for neuron_i = 1:size(mcc_map_info,1)
    
    neuralFilename = mcc_map_info.session{neuron_i};
    neuronLabel = mcc_map_info.unit{neuron_i};
    fprintf('Extracting: %s - %s ... [%i of %i]  \n',neuralFilename,neuronLabel,neuron_i,size(mcc_map_info,1))
   
    [h,p,~,stats] = ttest2(sdf_noncanc_all_saccade(neuron_i,timewins.error_comp+timewins.zero),...
        sdf_nostop_all_saccade(neuron_i,timewins.error_comp+timewins.zero));
    mean_noncanc = nanmean(sdf_noncanc_all_saccade(neuron_i,timewins.error_comp+timewins.zero));
    mean_nostop = nanmean(sdf_nostop_all_saccade(neuron_i,timewins.error_comp+timewins.zero));
    error_dir = mean_noncanc > mean_nostop;
    
    error_table(neuron_i,:) = table({neuralFilename},{neuronLabel},h,stats,error_dir,mean_noncanc,mean_nostop,...
        'VariableNames',{'neuralFilename','unit','h','stats','nc_dir','mean_nc','mean_ns'});
    
end

neuron_idx.error.nc_fac = find(error_table.h == 1 & error_table.nc_dir == 1);
neuron_idx.error.go_fac = find(error_table.h == 1 & error_table.nc_dir == 0);



%% Analyse: Clustering approach for errors
sdfWindow = timewins.sdf;

inputNeurons = []; inputNeurons = neuron_idx.error.nc_fac ;
inputSDF_target = {norm_sdf.target.noncanc(inputNeurons,:),norm_sdf.target.nostop(inputNeurons,:)};
inputSDF_error = {norm_sdf.saccade.noncanc(inputNeurons,:),norm_sdf.saccade.nostop(inputNeurons,:)};
inputSDF_tone = {norm_sdf.tone.noncanc(inputNeurons,:),norm_sdf.tone.nostop(inputNeurons,:)};

sdfTimes = {sdfWindow};
sdfEpoch = {[0:500]};

colorMapping = [1];

[sortIDs,idxDist, raw, respSumStruct, rawLink,myK] =...
    consensusCluster({norm_sdf.saccade.diff(inputNeurons,:)},sdfTimes,'-e',sdfEpoch,'-ei',colorMapping,'-er',sdfEpoch,'-c',0.5);

nClusters_manual = myK; clusterNeurons = [];
for cluster_i = 1:nClusters_manual
    clusterNeurons{cluster_i} = inputNeurons(sortIDs(:,nClusters_manual) == cluster_i );
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
colormap(viridis);
xlabel('Unit Number'); set(gca,'YAxisLocation','Left');
set(gca,'CLim',[-1 1])


% Generate a quick sdf of each cluster
for cluster_i = 1:nClusters_manual
    cluster_shift = 4*(cluster_i-1);
    figure('Renderer', 'painters', 'Position', [100 100 500 300]);hold on

    a = subplot(1,4,1); hold on
    plot(sdfTimes{1},nanmean(norm_sdf.target.noncanc(clusterNeurons{cluster_i},:),1), 'color', [colors.noncanc]);
    plot(sdfTimes{1},nanmean(norm_sdf.target.nostop(clusterNeurons{cluster_i},:),1), 'color', [colors.nostop]);
    vline([0], 'k--'); xlim([-500 250])

    b = subplot(1,4,[2 3 4]); hold on    
    plot(sdfTimes{1},nanmean(norm_sdf.saccade.noncanc(clusterNeurons{cluster_i},:),1), 'color', [colors.noncanc]);
    plot(sdfTimes{1},nanmean(norm_sdf.saccade.nostop(clusterNeurons{cluster_i},:),1), 'color', [colors.nostop]);
    

    vline([0 500 1000], 'k--'); xlim([-200 1500])

    linkaxes([a,b],'y');
    
    title(['Cluster ' int2str(cluster_i) ' - n: ' int2str(length(clusterNeurons{cluster_i}))])
    
end

%% Curation: Combine clusters
cluster_transient = [1,7,29,23,19,14];
cluster_sustained = [27,28,8,9,15,10,11];
cluster_noise = [2,3,4,5,6,12,13,16,17,18,20,21,24,25,26,22];

neuron_idx.transient = sort(cat(1,clusterNeurons{cluster_transient}));
neuron_idx.sustained = sort(cat(1,clusterNeurons{cluster_sustained}));
neuron_idx.noise = sort(cat(1,clusterNeurons{cluster_noise}));


































% 
% 
% 
% 
%% % Extract: Produce summary sheet figure for mean SDF
% for cluster_i = 1:nClusters_manual
% n_plot_x = 4; n_plot_y = 3; n_plot_sheet = n_plot_x*n_plot_y;
% n_batches = round(length(clusterNeurons{cluster_i})/n_plot_sheet,-1)+1;
% 
% neuron_list_i = 0;
% 
% for page_i = 1:n_batches
%     fig_out = figure('Renderer', 'painters', 'Position', [100 100 1200 800]);
% 
%     for plot_i = 1:n_plot_sheet
%         neuron_list_i = neuron_list_i+1;
%         try
%             neuron_i = clusterNeurons{cluster_i}(neuron_list_i);
%             subplot(n_plot_x, n_plot_y, plot_i); hold on
%            
%             
%             plot(timewins.sdf, norm_sdf.saccade.noncanc(neuron_i,:),'color',colors.noncanc)
%             plot(timewins.sdf, norm_sdf.saccade.nostop(neuron_i,:),'color',colors.nostop)
%             xlim([-200 600]); vline(0,'k--'); %, hline([-3 3],'r:'); hline([-6 6],'r:');
%             xlabel('Time from Target (ms)')
%             ylabel('Firing rate (z-score)')
%             title([mcc_map_info.session{neuron_i} ': ' mcc_map_info.unit{neuron_i}])
% 
% 
%         catch
%             continue
%         end
% 
%     end
% 
%     filename = fullfile(dirs.root,'results','sdf_overview_figs',['cluster_' int2str(cluster_i) '_errorsdf_saccade_overview_pg' int2str(page_i) '.pdf']);
%     set(fig_out,'PaperSize',[20 10]); %set the paper size to what you want
%     print(fig_out,filename,'-dpdf') % then print it
%     close(fig_out)
% end
% 
% end


%% Extract: Produce summary sheet figure for mean SDF
% 
% n_plot_x = 4; n_plot_y = 3; n_plot_sheet = n_plot_x*n_plot_y;
% n_batches = round(size(neuron_idx.error.nc_fac,1)/n_plot_sheet,-1)+1;
% 
% list_i = 0;
% 
% for page_i = 1:n_batches
%     fig_out = figure('Renderer', 'painters', 'Position', [100 100 1200 800]);
%     
%     for plot_i = 1:n_plot_sheet
%         list_i = list_i+1;
%         neuron_i = neuron_idx.error.nc_fac(list_i);
%         try
%             subplot(n_plot_x, n_plot_y, plot_i); hold on
%             plot(timewins.sdf, sdf_noncanc_all_saccade(neuron_i,:),'color',colors.noncanc)
%             plot(timewins.sdf, sdf_nostop_all_saccade(neuron_i,:),'color',colors.nostop)
%             xlim([-200 800]); vline(0,'k--'); %, hline([-3 3],'r:'); hline([-6 6],'r:'); 
%             xlabel('Time from Saccade (ms)')
%             ylabel('Firing rate (z-score)')
%             title([mcc_map_info.session{neuron_i} ': ' mcc_map_info.unit{neuron_i}])
%             
%             
%         catch
%             continue
%         end
%         
%     end
%     
%     filename = fullfile(dirs.root,'results','sdf_overview_figs',['errorsdf_saccade_overview_pg' int2str(page_i) '.pdf']);
%     set(fig_out,'PaperSize',[20 10]); %set the paper size to what you want
%     print(fig_out,filename,'-dpdf') % then print it
%     close(fig_out)
% end

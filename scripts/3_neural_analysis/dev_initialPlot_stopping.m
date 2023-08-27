
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
    ssd_act_canc = []; ssd_act_nostop = []; 
    % Latency match
    for ssd_i = 1:length(behavior(behaviorIdx,:).stopSignalBeh.inh_SSD)
        trl_canceled = []; trl_canceled = behavior(behaviorIdx,:).ttm.C.C{ssd_i};
        trl_nostop = []; trl_nostop = behavior(behaviorIdx,:).ttm.C.GO{ssd_i};
        
        if length(trl_canceled) < 10 | length(trl_nostop) < 10
            sdf_canceled_ssdx(ssd_i,:) = nan(1,length(-1000:2000));
            sdf_nostop_ssdx(ssd_i,:) = nan(1,length(-1000:2000));
            
            ssd_act_nostop(1,ssd_i) = nan;
            ssd_act_canc(1,ssd_i) = nan;
            
        else
            sdf_canceled_ssdx(ssd_i,:) = nanmean(data_in.SDF.stopSignal_artifical(trl_canceled,:));
            sdf_nostop_ssdx(ssd_i,:) = nanmean(data_in.SDF.stopSignal_artifical(trl_nostop,:));
            
            ssd_act_nostop(1,ssd_i) = nanmean(nanmean(data_in.SDF.stopSignal_artifical(trl_nostop,1000+[0:500])));
            ssd_act_canc(1,ssd_i) = nanmean(nanmean(data_in.SDF.stopSignal_artifical(trl_canceled,1000+[0:500])));        
        end
    end
    
    sdf_canceled_all_stopsignal(neuron_i,:) = nanmean(sdf_canceled_ssdx);
    sdf_nostop_all_stopsignal(neuron_i,:) = nanmean(sdf_nostop_ssdx);
    
    ssd_neurometric_nostop{neuron_i} = ssd_act_nostop;
    ssd_neurometric_canc{neuron_i} = ssd_act_canc;
    
end

%% Extract: Produce summary sheet figure for mean SDF

n_plot_x = 4; n_plot_y = 3; n_plot_sheet = n_plot_x*n_plot_y;
n_batches = round(size(mcc_map_info,1)/n_plot_sheet,-1)+1;

neuron_i = 0;
for page_i = 1:n_batches
    fig_out = figure('Renderer', 'painters', 'Position', [100 100 1200 800]);
    
    for plot_i = 1:n_plot_sheet
        neuron_i = neuron_i+1;
        
        neuralFilename = mcc_map_info.session{neuron_i};
        neuronLabel = mcc_map_info.unit{neuron_i};
        
        %... and find the corresponding behavior file index
        behFilename = data_findBehFile(neuralFilename);
        behaviorIdx = find(strcmp(dataFiles_beh,behFilename(1:end-4)));
        
        try
            subplot(n_plot_x, n_plot_y, plot_i); hold on
            plot(-1000:2000, sdf_canceled_all_stopsignal(neuron_i,:),'r')
            plot(-1000:2000, sdf_nostop_all_stopsignal(neuron_i,:),'k')
            xlim([-200 1000]); vline(0,'k-')
            xlabel('Time from stop-signal (ms)')
            ylabel('Firing rate (spks/sec)')
            title(['Neuron: ' int2str(neuron_i) ])
            vline(behavior.stopSignalBeh(behaviorIdx).ssrt.integrationWeighted,'k--')
        catch
            continue
        end
        
    end
    
    filename = fullfile(dirs.root,'results','sdf_overview_figs',['sdf_ssd_overview_pg' int2str(page_i) '.pdf']);
    set(fig_out,'PaperSize',[20 10]); %set the paper size to what you want
    print(fig_out,filename,'-dpdf') % then print it
    close(fig_out)
end



%% Extract: Produce summary sheet figure for neurometrics

n_plot_x = 4; n_plot_y = 3; n_plot_sheet = n_plot_x*n_plot_y;
n_batches = round(size(mcc_map_info,1)/n_plot_sheet,-1)+1;

neuron_i = 0;
for page_i = 1:n_batches
    fig_out = figure('Renderer', 'painters', 'Position', [100 100 1200 800]);
    
    for plot_i = 1:n_plot_sheet
        neuron_i = neuron_i+1;
        
        neuralFilename = mcc_map_info.session{neuron_i};
        neuronLabel = mcc_map_info.unit{neuron_i};
        
        %... and find the corresponding behavior file index
        behFilename = data_findBehFile(neuralFilename);
        behaviorIdx = find(strcmp(dataFiles_beh,behFilename(1:end-4)));
        
        try
            subplot(n_plot_x, n_plot_y, plot_i); hold on
            plot(behavior.stopSignalBeh(behaviorIdx).inh_SSD , ssd_neurometric_canc{neuron_i},'r','LineWidth',2)
            plot(behavior.stopSignalBeh(behaviorIdx).inh_SSD, ssd_neurometric_nostop{neuron_i},'k','LineWidth',2)
            
            scatter(behavior.stopSignalBeh(behaviorIdx).inh_SSD , ssd_neurometric_canc{neuron_i},'r','filled')
            scatter(behavior.stopSignalBeh(behaviorIdx).inh_SSD , ssd_neurometric_nostop{neuron_i},'k','filled')

            
            xlim([0 600])
            xlabel('Stop-signal delay (ms)')
            ylabel('Firing rate (spks/sec)')
            title(['Neuron: ' int2str(neuron_i) ])
        catch
            continue
        end
        
    end
    
    filename = fullfile(dirs.root,'results','sdf_overview_figs',['neurometric_ssd_overview_pg' int2str(page_i) '.pdf']);
    set(fig_out,'PaperSize',[20 10]); %set the paper size to what you want
    print(fig_out,filename,'-dpdf') % then print it
    close(fig_out)
end

% Difference
n_plot_x = 4; n_plot_y = 3; n_plot_sheet = n_plot_x*n_plot_y;
n_batches = round(size(mcc_map_info,1)/n_plot_sheet,-1)+1;

neuron_i = 0;
for page_i = 1:n_batches
    fig_out = figure('Renderer', 'painters', 'Position', [100 100 1200 800]);
    
    for plot_i = 1:n_plot_sheet
        neuron_i = neuron_i+1;
        
        neuralFilename = mcc_map_info.session{neuron_i};
        neuronLabel = mcc_map_info.unit{neuron_i};
        
        %... and find the corresponding behavior file index
        behFilename = data_findBehFile(neuralFilename);
        behaviorIdx = find(strcmp(dataFiles_beh,behFilename(1:end-4)));
        
        try
            subplot(n_plot_x, n_plot_y, plot_i); hold on
            plot(behavior.stopSignalBeh(behaviorIdx).inh_SSD , ssd_neurometric_canc{neuron_i}-ssd_neurometric_nostop{neuron_i},'b','LineWidth',2)
            scatter(behavior.stopSignalBeh(behaviorIdx).inh_SSD , ssd_neurometric_canc{neuron_i}-ssd_neurometric_nostop{neuron_i},'b','filled')

            
            xlim([0 600])
            xlabel('Stop-signal delay (ms)')
            ylabel('Firing rate (spks/sec)')
            title(['Neuron: ' int2str(neuron_i) ])
        catch
            continue
        end
        
    end
    
    filename = fullfile(dirs.root,'results','sdf_overview_figs',['neurometric_ssd_overview_diff_pg' int2str(page_i) '.pdf']);
    set(fig_out,'PaperSize',[20 10]); %set the paper size to what you want
    print(fig_out,filename,'-dpdf') % then print it
    close(fig_out)
end
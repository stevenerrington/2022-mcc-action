%% Place neurons into structure divided by type
neuron_index.ssd = find(mcc_analysis_table.glm_ssd_canc == 1);

input_neurons = [];
input_neurons = neuron_index.ssd;

n_plot_x = 4; n_plot_y = 3; n_plot_sheet = n_plot_x*n_plot_y;
n_batches = ceil(size(input_neurons,1)/n_plot_sheet);

%% Plot spike density function
neuron_i = 0;
for page_i = 1:n_batches
    fig_out = figure('Renderer', 'painters', 'Position', [100 100 1200 800]);
    
    for plot_i = 1:n_plot_sheet
        neuron_i = neuron_i+1;
        try
            neuron_j = input_neurons(neuron_i);
            
            neuralFilename = mcc_map_info.session{neuron_j};
            neuronLabel = mcc_map_info.unit{neuron_j};
            
            %... and find the corresponding behavior file index
            behFilename = data_findBehFile(neuralFilename);
            behaviorIdx = find(strcmp(dataFiles_beh,behFilename(1:end-4)));
            
            
            subplot(n_plot_x, n_plot_y, plot_i); hold on
            
            input_sdf_canc = [];
            input_sdf_canc = mcc_analysis_table.sdf_canceled{neuron_j}...
                (behavior(behaviorIdx,:).stopSignalBeh.midSSDidx+[-1, 0, +1],:);
            
            input_sdf_nostop = [];
            input_sdf_nostop = mcc_analysis_table.sdf_nostop{neuron_j}...
                (behavior(behaviorIdx,:).stopSignalBeh.midSSDidx+[-1, 0, +1],:);
            
            plot(-1000:2000, input_sdf_canc(1,:),'color',[colors.canceled 0.33])
            plot(-1000:2000, input_sdf_canc(2,:),'color',[colors.canceled 0.66])
            plot(-1000:2000, input_sdf_canc(3,:),'color',[colors.canceled 0.99])
            
%             plot(-1000:2000, input_sdf_nostop(1,:),'color',[colors.nostop 0.33])
%             plot(-1000:2000, input_sdf_nostop(2,:),'color',[colors.nostop 0.66])
%             plot(-1000:2000, input_sdf_nostop(3,:),'color',[colors.nostop 0.99])
            
            xlim([-250 750]); vline(0,'k-')
            xlabel('Time from stop-signal (ms)')
            ylabel('Firing rate (spks/sec)')
            title(['Neuron: ' int2str(neuron_j) ])
            vline(behavior.stopSignalBeh(behaviorIdx).ssrt.integrationWeighted,'k--')
        catch
            continue
        end
        
    end
    
    filename = fullfile(dirs.root,'results','glm_ssd',['glm_ssd_' trial_type_i{1} '_' int2str(page_i) '_sdf.pdf']);
    set(fig_out,'PaperSize',[20 10]); %set the paper size to what you want
    print(fig_out,filename,'-dpdf') % then print it
    close(fig_out)
end

%% Plot neurometric
neuron_i = 0;
for page_i = 1:n_batches
    fig_out = figure('Renderer', 'painters', 'Position', [100 100 1200 800]);
    
    for plot_i = 1:n_plot_sheet
        neuron_i = neuron_i+1;
        try
            neuron_j = input_neurons(neuron_i);
            
            neuralFilename = mcc_map_info.session{neuron_j};
            neuronLabel = mcc_map_info.unit{neuron_j};
            
            %... and find the corresponding behavior file index
            behFilename = data_findBehFile(neuralFilename);
            behaviorIdx = find(strcmp(dataFiles_beh,behFilename(1:end-4)));
            
            
            subplot(n_plot_x, n_plot_y, plot_i); hold on
            
            input_sdf_canc = [];
            input_sdf_canc = mcc_analysis_table.sdf_canceled{neuron_j};

            input_sdf_valid_idx = [];
            input_sdf_valid_idx = find(~isnan(input_sdf_canc(:,1)));
              
            a = []; b = []; c = [];
            for ssd_i = 1:length(input_sdf_valid_idx)
                a(ssd_i) = behavior(behaviorIdx,:).stopSignalBeh.inh_SSD(input_sdf_valid_idx(ssd_i));
                b(ssd_i) = behavior(behaviorIdx,:).stopSignalBeh.inh_pnc(input_sdf_valid_idx(ssd_i));
                c(ssd_i) = nanmean(input_sdf_canc(input_sdf_valid_idx(ssd_i),1000+[0:600]));
            end
            
            plot(a,b,'color','k','LineWidth',2)
            scatter(a,b,'filled','MarkerFaceColor','k')
            xlim([0 500]); ylim([0 1]); ylabel('p(respond | stop-signal)')
            
            yyaxis right
            
            plot(a,c,'color',colors.noncanc,'LineWidth',2)
            scatter(a,c,'filled','MarkerFaceColor',colors.noncanc)
            xlim([0 500]); ylabel('p(respond | stop-signal)')

            xlabel('Stop-signal delay (ms)')
            title(['Neuron: ' int2str(neuron_j) ])
        catch
            continue
        end
        
    end
    
    filename = fullfile(dirs.root,'results','glm_ssd',['glm_ssd_neurometric_' trial_type_i{1} '_' int2str(page_i) '_sdf.pdf']);
    set(fig_out,'PaperSize',[20 10]); %set the paper size to what you want
    print(fig_out,filename,'-dpdf') % then print it
    close(fig_out)
end

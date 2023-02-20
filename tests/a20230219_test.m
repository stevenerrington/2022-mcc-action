session_i = 15;

filename_neural = dajo_datamap_curated.dataFilename{session_i};
filename_beh = data_findBehFile(filename_neural);
index_beh = find(strcmp(behavior.sessionName,filename_beh));

clear spikes
load(fullfile(dirs.data,[filename_neural '-spk.mat']))

[sdf, spk, ~] = spk_alignTrials(behavior.trialEventTimes{index_beh}(:,[2,3,5,6,9]), spikes.time, [-1000 2500]);


%% Generate plot using gramm


n_neurons = size(dajo_datamap_curated.spkInfo(session_i,1).unitInfo,1);
xlim_input = [-200 1500];
trial_epoch_list = {'target','saccade','tone'};

for neuron_i = 1:n_neurons
    
    neuron_label = dajo_datamap_curated.spkInfo(session_i,1).unitInfo.unitDSP{neuron_i};
    
    clear figure_plot trials_in
    trials_in = behavior.ttx(index_beh).nostop.all.all;
    
    for epoch_i = 1:length(trial_epoch_list)
        trial_epoch = trial_epoch_list{epoch_i};

        % Raster plot
        figure_plot(1,epoch_i)=gramm('x',spk.(neuron_label).(trial_epoch)(trials_in));
        figure_plot(1,epoch_i).geom_raster();
        figure_plot(1,epoch_i).axe_property('XLim',xlim_input);
        figure_plot(1,epoch_i).set_names('x',['Time from ' trial_epoch ' (ms)'],'y','Trial');
        figure_plot(1,epoch_i).geom_vline('xintercept',0,'style','k-');
        
        % Spike density function
        figure_plot(2,epoch_i)=gramm('x',[-1000:2500],'y',num2cell(sdf.(neuron_label).(trial_epoch)(trials_in,:),2));
        figure_plot(2,epoch_i).stat_summary();
        figure_plot(2,epoch_i).axe_property('XLim',xlim_input);
        figure_plot(2,epoch_i).set_names('x',['Time from ' trial_epoch ' (ms)'],'y','Firing rate (spk/sec)');
        figure_plot(2,epoch_i).geom_vline('xintercept',0,'style','k-');
        
    end
    
    
    
    figure_plot_out = figure('Renderer', 'painters', 'Position', [100 100 1000 700]);
    figure_plot.draw();
    
end

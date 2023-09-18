clear plot_*

%% Extract activity within a window, following stop-signal delay

% For each modulated neuron
for neuron_idx = 1:size(mcc_analysis_table,1)
    
    % Get the corresponding file name 
    neuralFilename = mcc_analysis_table.session{neuron_idx};
    neuronLabel = mcc_analysis_table.unit{neuron_idx};
    fprintf('Extracting: %s - %s ...   \n',neuralFilename,neuronLabel)
    
    %... and find the corresponding behavior file index
    behFilename = data_findBehFile(neuralFilename);
    behaviorIdx = find(strcmp(dataFiles_beh,behFilename(1:end-4)));
    
    % Find the middle most SSD [-1 0 1 2];
    ssd_idx = []; ssd_idx = behavior(behaviorIdx,:).stopSignalBeh.midSSDidx + [-1 0 1 2];
    
    % Save the output for each neuron
    plot_ssd_out(neuron_idx,:) = behavior(behaviorIdx,:).stopSignalBeh.inh_SSD(ssd_idx);
    plot_ssdidx_out(neuron_idx,:) = [-1 0 1 2];
    plot_pnc_out(neuron_idx,:) = behavior(behaviorIdx,:).stopSignalBeh.inh_pnc(ssd_idx);
    
    average_fr = [];
    average_fr = nanmean(mcc_analysis_table.sdf_canceled{neuron_idx}...
        (ssd_idx,1000+[0:600]),2);
    
    plot_fr_out(neuron_idx,:) = abs(average_fr)./max(abs(average_fr))';
    
end

%% Plot average windowed activity as a function of stop-signal delay and p(responding)
clear figure_plot

% Loop through all the neurons, concatenate their indices, and get the
% corresponding cluster label for plotting
neurons_in = [];
cluster_label = [];
for cluster_i = 1:length(cluster_merge_idx)
    neurons_in = [neurons_in; cluster_neuron_id_table{cluster_i}];
    cluster_label = [cluster_label; repmat({[int2str(cluster_i) '_cluster']},length(cluster_neuron_id_table{cluster_i}),1)];
end

% Produce the figure in gramm
figure_plot(1,1)=gramm('x',[plot_ssdidx_out(neurons_in,:); plot_ssdidx_out(neurons_in,:)],...
    'y',[plot_pnc_out(neurons_in,:); plot_fr_out(neurons_in,:)],...
    'color',[repmat({'0_pnc'},length(neurons_in),1); cluster_label]);
figure_plot(1,1).stat_summary('geom',{'line','point','errorbar'});
figure_plot(1,1).axe_property('YLim',[0 1.0]);
figure_plot(1,1).set_names('y','');

% Draw the figure
figure('Renderer', 'painters', 'Position', [100 100 400 400]);
figure_plot.draw

%% Plot the spike density function for each cluster, split by ranked stop-signal delay
clear sdf_plot

% Define a colormap that is sequential
colormap_canc = colormap(cbrewer('seq', 'Reds', 8));
colormap_canc = colormap_canc(3:6,:);

% For each cluster
for cluster_i = 1:length(cluster_merge_idx)
    
    % Define the neurons to look at 
    input_neurons = []; input_neurons = cluster_neuron_id_table{cluster_i};
    
    % Initialize the arrays (as we will concatenate into this)
    canceled_sdf = [];
    canceled_sdf_label = [];
    
    % For each of our neurons
    for neuron_i = 1:length(input_neurons)
        % Get the actual index in the main table
        neuron_j = input_neurons(neuron_i);
        
        % Get the filename and information
        neuralFilename = mcc_analysis_table.session{neuron_j};
        neuronLabel = mcc_analysis_table.unit{neuron_j};
        fprintf('Extracting: %s - %s ...   \n',neuralFilename,neuronLabel)
        
        %... and find the corresponding behavior file index
        behFilename = data_findBehFile(neuralFilename);
        behaviorIdx = find(strcmp(dataFiles_beh,behFilename(1:end-4)));
        
        % Get the index of the middlemost SSD [-1 0 1 2]
        ssd_idx = []; ssd_idx = behavior(behaviorIdx,:).stopSignalBeh.midSSDidx + [-1 0 1 2];
         
        average_fr = [];
        average_fr = mcc_analysis_table.sdf_canceled{neuron_j}(ssd_idx,:);
        
        norm_fr_max = []; norm_fr_min = []; norm_fr = [];
        norm_fr_max = max(average_fr(:,1000+[0:600]),[],2);
        norm_fr_min = min(average_fr(:,1000+[0:600]),[],2);
        
        norm_fr = [norm_fr_max; norm_fr_min];
        index = []; [~, index] = max(abs(norm_fr));

        % Get the SDF for the given neuron in this index
        canceled_sdf = [canceled_sdf;...
            average_fr./abs(norm_fr(index))];
        
        % And produce the relevant label for plotting
        canceled_sdf_label = [canceled_sdf_label; {'1','2','3','4'}'];
        
    end
    
    % For the given cluster, produce the SDF figure
    sdf_plot(1,cluster_i)=gramm('x',[-1000:2000],...
        'y',canceled_sdf,...
        'color',canceled_sdf_label);
    
    % As a summary SDF (-/+ SEM)
    sdf_plot(1,cluster_i).stat_summary();
    sdf_plot(1,cluster_i).axe_property('XLim',[-200 600],'YLim',[-0.5 0.5]);
    sdf_plot(1,cluster_i).set_names('y','');
    sdf_plot(1,cluster_i).no_legend();
    sdf_plot(1,cluster_i).set_color_options('map',colormap_canc);
    
end

% Draw the SDF figure
figure('Renderer', 'painters', 'Position', [100 100 1000 250]);
sdf_plot.draw

%% Plot the proportion of neurons that had significant periods of modulation for SSD
%  and value, based on the GLM

% Open the figure window
figure('Renderer', 'painters', 'Position', [100 100 1000 500]);

% For each cluster
for cluster_i = 1:5
    % Define the input window
    input_neurons = []; input_neurons = cluster_neuron_id_table{cluster_i};
    
    % Find the number of neurons that had a sig GLM hit for SSD
    n_sig_ssd = sum(mcc_analysis_table.glm_ssd_canc(input_neurons));
    n_nonsig_ssd = length(input_neurons)-n_sig_ssd;
    
    % Find the number of neurons that had a sig GLM hit for Value
    n_sig_value = sum(mcc_analysis_table.glm_value_canc(input_neurons));
    n_nonsig_value = length(input_neurons)-n_sig_value;
    
    % Plot SSD hits as a donut plot
    subplot(2,5,cluster_i);
    donut([n_sig_ssd, n_nonsig_ssd]);
    
    % Plot value hits as a donut plot
    subplot(2,5,cluster_i+5);
    donut([n_sig_value, n_nonsig_value]);
    
end



